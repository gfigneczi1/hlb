function [segments] = segmentor_mapValidation(input_folder)
% Segmentation rules
% 1. Consider only RTK_Fixed parts of the map, but with debouncing!
% Therefore if the RTK_Fixed disappears to a specific time, but returns
% within an interval, the segment shall not be cut!
    zalaData = 0;

    matFiles = dir(fullfile(input_folder,'/*.mat'));
    matFilesTable = struct2table(matFiles);
    clear refMatFiles;
    sortedMatFilesTable = sortrows(matFilesTable, 'date');
    matFiles = table2struct(sortedMatFilesTable);
        
    data = cell(1,length(matFiles));
    testFileNumber = 1;
    mandatorySignals = ["GPS_status", "VelocityX_ESP", "GPS_time", "LongPos_abs", "LatPos_abs"];
    testSignals = ["c01_left","c01_right","c1","c2"];
    
    for fileID = 1:length(matFiles)
        rawData = load(fullfile(matFiles(fileID).folder,matFiles(fileID).name));
        if (isfield(rawData, 'rawData'))
            rawData = rawData.rawData;
        elseif (isfield(rawData, 'segment'))
            rawData = rawData.segment;
        end
        if (isfield(rawData, "LongPos_Abs"))
            rawData.LongPos_abs = rawData.LongPos_Abs;            
            rawData=rmfield(rawData, 'LongPos_Abs');
        end
        if (isfield(rawData, "LatPos_Abs"))
            rawData.LatPos_abs = rawData.LatPos_Abs;
            rawData=rmfield(rawData, 'LatPos_Abs');
        end
        [ mandatoryStatus, testStatus ] = signalChecker( rawData, mandatorySignals, testSignals );
        if (zalaData && mandatoryStatus == 0)
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
        end
        [ mandatoryStatus, testStatus ] = signalChecker( rawData, mandatorySignals, testSignals );
            
        disp("Result of signal checker:");
        disp(strcat("mandatory signals:",num2str(mandatoryStatus)));
        disp(strcat("test signals:",num2str(testStatus)));
        if (testStatus == 1)
            % VALIDATOR FILE
            [rawData.X_abs,rawData.Y_abs] = deg2utm(rawData.LatPos_abs, rawData.LongPos_abs);
            rawData.X_abs = rawData.X_abs';
            rawData.Y_abs = rawData.Y_abs';
            rawData = redundantPointsFilter(rawData, "Time");
            if (contains(matFiles(fileID).name, "south"))
                rawData.theta_calc = theta_recalc(rawData,1);
            else
                rawData.theta_calc = theta_recalc(rawData,0);
            end
            rawData.invalidVx = rawData.VelocityX_ESP < 10;
            rawData.validc0 = abs(0.5*(rawData.c01_left+rawData.c01_right)) < 1.5; 
            [rawData.localization, rawData.roadType] = localization(rawData, input_folder);
            rawData = debouncer(rawData, "GPS_status", 8, 20); % last arg is debounce time in [s]
            rawData = structure_filter(rawData, "GPS_status", 4, "GT");
            %rawData = structure_filter(rawData, "lanevalidity", 1);
            % 1: rural road, 2: highway
            rawData = structure_filter(rawData, "roadType", 1, "EQ");
            rawData = structure_filter(rawData, "invalidVx",0, "EQ");
            % only necessary if opposite direction traffic is supposed
            rawData = structure_filter(rawData, "validc0",1, "EQ");
        elseif (mandatoryStatus == 1)
            % only mandatory fields, no test signals
            [rawData.X_abs,rawData.Y_abs] = deg2utm(rawData.LatPos_abs, rawData.LongPos_abs);
            rawData.X_abs = rawData.X_abs';
            rawData.Y_abs = rawData.Y_abs';
            rawData = redundantPointsFilter(rawData, "Time");
            if (contains(matFiles(fileID).name, "south"))
                rawData.theta_calc = theta_recalc(rawData,1);
            else
                rawData.theta_calc = theta_recalc(rawData,0);
            end
            rawData.invalidVx = rawData.VelocityX_ESP < 10;
            [rawData.localization, rawData.roadType] = localization(rawData, input_folder);
            rawData = debouncer(rawData, "GPS_status", 4, 20); % last arg is debounce time in [s]
            rawData = structure_filter(rawData, "GPS_status", 4, "GT");
            % 1: rural road, 2: highway
            rawData = structure_filter(rawData, "roadType", 1, "EQ");
            rawData = structure_filter(rawData, "invalidVx",0, "EQ");
        end
        signals = fieldnames(rawData);
        n = length(signals);
        m = size(eval(['rawData.' signals{1}]),2);
        matFileData = reshape(table2array(struct2table(rawData)),[m n]);
        data{fileID} = matFileData;
        data_table = array2table(vertcat(data{fileID}), 'VariableNames',signals);
        if (testStatus == 1)
            if (exist('segments')) 
                if (~isfield(segments,'validator'))
                    segments.validator.segment = data_table;
                    segments.validator.name = matFiles(fileID).name(1:min(length(matFiles(fileID).name),50));
                else
                    disp('Multiple validator files are found. Only the first one is added!');
                end
            else
                segments.validator.segment = data_table;
                segments.validator.name = matFiles(fileID).name(1:min(length(matFiles(fileID).name),50));
            end
        else
            segments.segments(testFileNumber).segment = data_table;
            segments.segments(testFileNumber).name = matFiles(fileID).name(1:min(length(matFiles(fileID).name),100));
            segments.segments(testFileNumber).signalStatus = [ mandatoryStatus, testStatus ];
            testFileNumber = testFileNumber + 1;
        end
        clear rawData;
        clear matFileData;
    end

    clear data data_table;
    
end