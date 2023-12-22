function [segments] = segmentor_OL_LDM_resim(input_folder)
% Segmentation rules
% 1. GPS_Time based redundant points filter (dt == 50ms)
% 2. lane validity for proper corridor and Lane Change detection (robust
% detection)
    matFiles = dir(fullfile(input_folder,'/*.mat'));
    MatFilesTable = struct2table(matFiles);
    clear matFiles;
    sortedMatFilesTable = sortrows(MatFilesTable, 'date');
    matFiles = table2struct(sortedMatFilesTable);
    
    data = cell(1,length(matFiles));
    for fileID = 1:length(matFiles)
        rawData = load(fullfile(matFiles(fileID).folder,matFiles(fileID).name));
%         rawData = rmfield(rawData, "q_KEY");
%         rawData = rmfield(rawData, "q_UNIT");
%         rawData = rmfield(rawData, "q_SCMIN");
%         rawData = rmfield(rawData, "q_SCMAX");
%         rawData = rmfield(rawData, "q_TC");
%         rawData = redundantPointsFilter(rawData, "Time");
%         rawData.X_abs = rawData.LongPos_abs * 40075000 .* cos(rawData.LatPos_abs*pi()/180) / 360;
%         rawData.Y_abs = rawData.LatPos_abs * 111.32*1000;
%         rawData.theta_calc = theta_recalc(rawData);
        % Restructuring
        signals = fieldnames(rawData);
        n = length(signals);
        m = size(eval(['rawData.' signals{1}]),2);
        matFileData = reshape(table2array(struct2table(rawData)),[m n]);
        data{fileID} = matFileData;
        data_subtable = array2table(matFileData, 'VariableNames',signals);
        segments.segments(fileID).segment = data_subtable;
        clear rawData;
        clear matFileData;
        clear data_subtable;
    end
    clear data data_table;    
end