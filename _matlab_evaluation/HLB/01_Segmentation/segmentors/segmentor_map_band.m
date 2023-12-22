function [segments] = segmentor_map_band(input_folder)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    matFiles = dir(fullfile(input_folder,'/*.mat'));
    MatFilesTable = struct2table(matFiles);
    clear matFiles;
    sortedMatFilesTable = sortrows(MatFilesTable, 'date');
    matFiles = table2struct(sortedMatFilesTable);
    testSignals = ["c01_left", "c01_right","c1","c2", "Right_Index", "Left_Index"]; % signals only available in a test vehicle
    % note: test signals can be stubbed by the converter by e.g. only
    % zeros, then the signal checker here will fail to detect missing
    % signals
    mandatorySignals = ["VelocityX_ESP", "GPS_status", "LongPos_abs", "LatPos_abs", "GPS_time"];
    
    data = cell(1,length(matFiles));
    
    for fileID = 1:length(matFiles)
        rawData = load(fullfile(matFiles(fileID).folder,matFiles(fileID).name));
        if (isfield(rawData,'segment'))
            rawData = rawData.segment;
        elseif (isfield(rawData, 'rawData'))
            rawData = rawData.rawData;
        end
        [ ~, testStatus ] = signalChecker( rawData, mandatorySignals, testSignals );
        if (testStatus == 1)
            [rawData.lanevalidity, rawData.LCL, rawData.LCR] = lanevaliditycheck(rawData);
            rawData.theta_calc = theta_recalc(rawData);
            %rawData = structure_filter(rawData, "GPS_status", 8, "EQ");
            signals = fieldnames(rawData);
            n = length(signals);
            m = size(eval(['rawData.' signals{1}]),2);
            matFileData = reshape(table2array(struct2table(rawData)),[m n]);
            data{1} = matFileData;
            clear rawData;
            clear matFileData;
            data_table = array2table(vertcat(data{:}), 'VariableNames', signals);
            segments = data_table;
        end
    end
    clear data data_table;
end

