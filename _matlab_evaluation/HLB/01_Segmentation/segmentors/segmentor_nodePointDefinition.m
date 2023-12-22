function [segments, config] = segmentor_nodePointDefinition(input_folder, config)
% Segmentation rules
% 1. GPS_Time based redundant points filter (dt == 50ms)
% 2. lane validity for proper corridor and Lane Change detection (robust
% detection)
    roadType = 1;

    matFiles = dir(fullfile(input_folder,'/*.mat'));
    MatFilesTable = struct2table(matFiles);
    clear matFiles;
    sortedMatFilesTable = sortrows(MatFilesTable, 'date');
    matFiles = table2struct(sortedMatFilesTable);
    
    data = cell(1,length(matFiles));
    for fileID = 1:length(matFiles)
        rawData = load(fullfile(matFiles(fileID).folder,matFiles(fileID).name));
        
        rawData = redundantPointsFilter(rawData, "Time");
        % Generic validation and cut
        [rawData.lanevalidity, rawData.LCL, rawData.LCR] = lanevaliditycheck(rawData);
        
        rawData.theta_calc = theta_recalc(rawData);
        [rawData.localization, rawData.roadType] = localization(rawData);
        
        rawData = debouncer(rawData, "GPS_status", 8, 20); % last arg is debounce time in [s]
        rawData = structure_filter(rawData, "GPS_status", 8);
        rawData = structure_filter(rawData, "lanevalidity", 1);
        rawData = structure_filter(rawData, "roadType", roadType);
        
        % Restructuring
        signals = fieldnames(rawData);
        n = length(signals);
        m = size(eval(['rawData.' signals{1}]),2);
        matFileData = reshape(table2array(struct2table(rawData)),[m n]);
        data{fileID} = matFileData;
        data_subtable = array2table(matFileData, 'VariableNames',signals);
        segments.segments(fileID).segment = data_subtable;
        segments.segments(fileID).meta.roadType = roadType;
        segments.segments(fileID).meta.file = matFiles(fileID).name;
        clear rawData;
        clear matFileData;
        clear data_subtable;
    end
    data_table = array2table(vertcat(data{:}), 'VariableNames',signals);
    segments.uncut = data_table;
    clear data data_table;
    
end