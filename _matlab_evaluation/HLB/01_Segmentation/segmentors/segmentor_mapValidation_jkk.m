function [segments] = segmentor_mapValidation_jkk(input_folder)
% Segmentation rules
% 1. Consider only RTK_Fixed parts of the map, but with debouncing!
% Therefore if the RTK_Fixed disappears to a specific time, but returns
% within an interval, the segment shall not be cut!

    matFiles = dir(fullfile(input_folder,'/zala*.mat'));
    matFilesTable = struct2table(matFiles);
    clear refMatFiles;
    sortedMatFilesTable = sortrows(matFilesTable, 'date');
    matFiles = table2struct(sortedMatFilesTable);
    
    mapFiles = dir(fullfile(input_folder,'/map*.mat'));
    if (length(mapFiles) ==0 )
        disp('Insufficient map files are found');
        segments.maps = struct();
    else
        % processing the map file
        for i=1:length(mapFiles)
            map = load(fullfile(mapFiles(i).folder, mapFiles(i).name));
            map = map.map;
            if contains(mapFiles(i).name,'north')
                if (contains(mapFiles(i).name, 'right'))
                    segments.maps.north_right = map;
                else
                    segments.maps.north_left = map;
                end
            else
                if (contains(mapFiles(i).name, 'right'))
                    segments.maps.south_right = map;
                else
                    segments.maps.south_left = map;
                end
            end
        end
    end
        
    data = cell(1,length(matFiles));
    testFileNumber = 1;
    mandatorySignals = ["GPS_status", "GPS_time", "LongPos_abs", "LatPos_abs"];
    testSignals = [];
    
    for fileID = 1:length(matFiles)
        rawData = load(fullfile(matFiles(fileID).folder,matFiles(fileID).name));
        if (isfield(rawData, 'rawData'))
            rawData = rawData.rawData;
        elseif (isfield(rawData, 'segment'))
            rawData = rawData.segment;
        end

        [ mandatoryStatus, testStatus ] = signalChecker( rawData, mandatorySignals, testSignals );
        disp("Result of signal checker:");
        disp(strcat("mandatory signals:",num2str(mandatoryStatus)));
        if (mandatoryStatus == 0)
            for i=1:size(rawData.X_abs_GNSS,2)
                UTC(i,1:4) = '33 T';
            end
            [rawData.LatPos_abs,rawData.LongPos_abs] = utm2deg(rawData.X_abs_GNSS',rawData.Y_abs_GNSS',UTC);
            clear UTC
            rawData.LatPos_abs = rawData.LatPos_abs';
            rawData.LongPos_abs = rawData.LongPos_abs';
            rawData.X_abs = rawData.X_abs_GNSS;
            rawData.Y_abs = rawData.Y_abs_GNSS;
            tStartIdx = find(rawData.GPS_time>0,1);
            rawData.GPS_time(tStartIdx(1,1):end) = rawData.GPS_time(tStartIdx(1,1):end) + rawData.q_T0(tStartIdx(1,1):end);
            rawData.Left_Index = zeros(size(rawData.X_abs));
            rawData.Right_Index = zeros(size(rawData.X_abs));
            rawData.yawRate_ESP = zeros(size(rawData.X_abs));
            
            rawData = redundantPointsFilter(rawData, "Time");
            rawData.theta_calc = theta_recalc(rawData);
            rawData.invalidVx = rawData.VelocityX_ESP < 10;
            
            [rawData.localization, rawData.roadType] = localization(rawData, input_folder);
            rawData = debouncer(rawData, "GPS_status", 4, 20); % last arg is debounce time in [s]
            rawData = structure_filter(rawData, "GPS_status", 4);
            rawData = structure_filter(rawData, "invalidVx",0);
            rawData = structure_filter(rawData, "localization",1);
            
            signals = fieldnames(rawData);
            n = length(signals);
            m = size(eval(['rawData.' signals{1}]),2);
            matFileData = reshape(table2array(struct2table(rawData)),[m n]);
            data{fileID} = matFileData;
            data_table = array2table(vertcat(data{fileID}), 'VariableNames',signals);
            
            segments.segments(testFileNumber).segment = data_table;
            segments.segments(testFileNumber).name = matFiles(fileID).name(1:min(length(matFiles(fileID).name),100));
            segments.segments(testFileNumber).signalStatus = [ mandatoryStatus, testStatus ];
            if (contains(matFiles(fileID).name,"undisturbed", 'IgnoreCase',true))
                if (contains(matFiles(fileID).name,"north", 'IgnoreCase',true))
                    segments.segments(testFileNumber).lane = "north_right";
                else
                    segments.segments(testFileNumber).lane = "south_right";
                end
            end                
            testFileNumber = testFileNumber + 1;
        end
        clear rawData;
        clear matFileData;
    end

    clear data data_table;
    
end