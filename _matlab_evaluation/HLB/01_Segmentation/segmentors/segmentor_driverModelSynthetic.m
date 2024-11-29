function [segments] = segmentor_driverModelSynthetic(input_folder)
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
    
    segmentID = 1;
    data = cell(1,length(matFiles));
    for fileID = 1:length(matFiles)
        load(fullfile(matFiles(fileID).folder,matFiles(fileID).name));
        
        % transforming to struct
        rawData.X_abs = data.egoPosition(:,1);
        rawData.Y_abs = data.egoPosition(:,2);
        rawData.LaneEdgePositionLeft = zeros(length(rawData.X_abs),1);
        rawData.LaneEdgePositionRight = zeros(length(rawData.X_abs),1);
        rawData.LaneOrientation = data.laneOrientation;
        rawData.LaneCurvature = data.input(:,6);
        rawData.c3 = data.input(:,7);
        rawData.VelocityX = data.input(:,3);
        rawData.AccelerationX = data.input(:,4);
        rawData.YawRate = data.input(:,5);
        rawData.OncomingTrafficType = data.input(:,1);
        rawData.FrontTrafficType = data.input(:,2);
        rawData.theta_calc = diff(rawData.Y_abs)./diff(rawData.X_abs);
        
        if (size(rawData.theta_calc,1) == 1)
            rawData.theta_calc = [rawData.theta_calc(1) rawData.theta_calc];
        else
            rawData.theta_calc = [rawData.theta_calc(1); rawData.theta_calc];
        end
        rawData.theta_calc = atan(rawData.theta_calc);
        if (mean(diff(rawData.X_abs))<0)
            rawData.theta_calc = rawData.theta_calc+pi();
        end
        rawData.Relative_time = data.relativeTime;

        clear data

        signals = fieldnames(rawData);
        n = length(signals);
        m = length(eval(['rawData.' signals{1}]));
        matFileData = reshape(table2array(struct2table(rawData)),[m n]);
        data{segmentID} = matFileData;
        data_subtable = array2table(matFileData, 'VariableNames',signals);
        segments.segments(segmentID).segment = data_subtable;
        segments.segments(segmentID).name = strcat(matFiles(fileID).name(1:min(length(matFiles(fileID).name),100)));
        clear rawData;
        clear matFileData;
        clear data_subtable;
        clear cutInfo;
        segmentID = segmentID + 1;
    end
    % This is commented out here: the uncut data contains all signals of
    % all matfiles, if matfiles are different, this causes error.
    clear data data_table;
    
end