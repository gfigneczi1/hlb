function [segments] = segmentor_driverModel(input_folder)
% Segmentation rules
% 1. GPS_Time based redundant points filter (dt == 50ms)
% 2. lane validity for proper corridor and Lane Change detection (robust
% detection)
% 3. Localization: consider only relevant positioned parts of the
% measurements

    matFiles = dir(fullfile(input_folder,'/*.mat'));
    MatFilesTable = struct2table(matFiles);
    clear matFiles;
    sortedMatFilesTable = sortrows(MatFilesTable, 'date');
    matFiles = table2struct(sortedMatFilesTable);
    testSignals = ["LaneEdgePositionLeft", ...
        "LaneEdgePositionRight", ...
        "LaneOrientation", ...
        "LaneCurvature"]; % signals only available in a test vehicle
    % note: test signals can be stubbed by the converter by e.g. only
    % zeros, then the signal checker here will fail to detect missing
    % signals
    mandatorySignals = ["VelocityX", "GPS_status", "LongPos_abs", "LatPos_abs", "GPS_time"];
    roadTypeFilter = 1;
      
    segmentID = 1;
    data = cell(1,length(matFiles));
    for fileID = 1:length(matFiles)
        lapnumber = 1;
        rawData = load(fullfile(matFiles(fileID).folder,matFiles(fileID).name));
        if (isfield(rawData,'segment'))
            rawData = rawData.segment;
        elseif (isfield(rawData, 'rawData'))
            rawData = rawData.rawData;
        end
        %[rawDatas, lapnumber] = circle_cutter(rawData);
        for i=1:lapnumber
            %rawData = rawDatas{i};
            [ mandatoryStatus, testStatus ] = signalChecker( rawData, mandatorySignals, testSignals );
            disp("Result of signal checker:");
            disp(strcat("mandatory signals:",num2str(mandatoryStatus)));
            disp(strcat("test signals:",num2str(testStatus)));
            if (mandatoryStatus && testStatus)
                rawData = redundantPointsFilter(rawData, "Time");
                % Generic validation and cut
                [rawData.localization, rawData.roadType] = localization(rawData, input_folder);        
                rawData = debouncer(rawData, "GPS_status", 8, 20); % last arg is debounce time in [s]
                rawData.invalidC0 = ((abs(0.5*(rawData.LaneEdgePositionLeft+rawData.LaneEdgePositionRight)) > 0.8)+(isnan(rawData.LaneEdgePositionLeft))+(rawData.LaneEdgePositionLeft==0)+(rawData.LaneEdgePositionLeft>3.5)+(rawData.LaneEdgePositionRight<-3.5)) > 0;
                rawData.invalidC1 = abs(rawData.LaneOrientation) > 0.05; 
                rawData.invalidSpeed = rawData.VelocityX < 5;
                rawData = extend(rawData,"invalidC0",1,-3);
                rawData = extend(rawData,"invalidC0",1,3);
                rawData = extend(rawData,"invalidC1",1,-3);
                rawData = extend(rawData,"invalidC1",1,3);
                
                rawData.isStraight = abs(movmean(rawData.LaneCurvature,100)) <= 2.5e-4;
                                
                rawData.X_abs = rawData.LongPos_abs * 40075000 .* cos(rawData.LatPos_abs*pi()/180) / 360;
                rawData.Y_abs = rawData.LatPos_abs * 111.32*1000;
                if (contains(matFiles(fileID).name, "62B") || contains(matFiles(fileID).name, "south"))
                    rawData.theta_calc = theta_recalc(rawData, 1);
                elseif (contains(matFiles(fileID).name, "M7A") && roadTypeFilter==2)
                    rawData.theta_calc = theta_recalc(rawData, 1);
                else
                    rawData.theta_calc = theta_recalc(rawData, 0);
                end
                
                uncut = rawData;
                cutInfo.roadType = rawData.roadType ~=1;
                cutInfo.invalidSpeed = rawData.invalidSpeed == 1;
                cutInfo.invalidC0 = rawData.invalidC0 == 1;
                cutInfo.GPS_status = rawData.GPS_status < 8;
                
                rawData = structure_filter(rawData, "GPS_status", 4, "GT");
                %rawData = structure_filter(rawData, "localization", 1, "EQ");
                % 1: country road; 2: highway
                %rawData = structure_filter(rawData, "roadType", roadTypeFilter, "EQ");
                %rawData = structure_filter(rawData, "isStraight", 1, "EQ");
                rawData = structure_filter(rawData, "invalidSpeed", 0, "EQ");
                %rawData = structure_filter(rawData, "invalidC0", 0, "EQ");
           
                %rawData = structure_filter(rawData, "invalidC1", 0, "EQ");
                if (isfield(rawData,'control_active'))
                    %rawData = structure_filter(rawData, "control_active", 3, "EQ");
                end
%                 if (~isfield(rawData,'vehicleOppositeSideDisturbing'))
%                     rawData.vehicleOppositeSideDisturbing = rawData.VelocityX_ESP * 0;
%                     rawData = extend(rawData, "vehicleOppositeSideDisturbing", -3);
%                     rawData = extend(rawData, "vehicleOppositeSideDisturbing", 1);
%                 end
                
            end
            if (testStatus)
                % Restructuring
                signals = fieldnames(rawData);
                n = length(signals);
                m = size(eval(['rawData.' signals{1}]),2);
                matFileData = reshape(table2array(struct2table(rawData)),[m n]);
                data{segmentID} = matFileData;
                data_subtable = array2table(matFileData, 'VariableNames',signals);
                segments.segments(segmentID).cutInfo = cutInfo;
                segments.segments(segmentID).segment = data_subtable;
                segments.segments(segmentID).name = strcat(matFiles(fileID).name(1:min(length(matFiles(fileID).name),100)),'_lap',num2str(i));
                segments.segments(segmentID).signalStatus = [mandatoryStatus testStatus];
                segments.segments(segmentID).uncut = uncut;
                clear rawData;
                clear matFileData;
                clear data_subtable;
                clear cutInfo;
                segmentID = segmentID + 1;
            end
        end
        if (mandatoryStatus && ~isfield(segments, 'uncut'))
            data_table = array2table(vertcat(data{:}), 'VariableNames',signals);
            segments.uncut = data_table;
        end
    end
    % This is commented out here: the uncut data contains all signals of
    % all matfiles, if matfiles are different, this causes error.
    clear data data_table;
    
end