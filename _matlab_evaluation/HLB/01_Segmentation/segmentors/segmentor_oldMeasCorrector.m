function [segments] = segmentor_oldMeasCorrector(input_folder, config)
% Segmentation rules
% 1. Consider only RTK_Fixed parts of the map, but with debouncing!
% Therefore if the RTK_Fixed disappears to a specific time, but returns
% within an interval, the segment shall not be cut!

    matFiles = dir(fullfile(input_folder,'/*.mat'));
    matFilesTable = struct2table(matFiles);
    clear refMatFiles;
    sortedMatFilesTable = sortrows(matFilesTable, 'date');
    matFiles = table2struct(sortedMatFilesTable);
        
    data = cell(1,length(matFiles));
    
    segments = [];
    
    for fileID = 1:length(matFiles)
        rawData = load(fullfile(matFiles(fileID).folder,matFiles(fileID).name));
        % calculating intermediate variables of absolute coordinates
        X_abs = rawData.LongPos_abs * 40075000 .* cos(rawData.LatPos_abs*pi()/180) / 360;
        Y_abs = rawData.LatPos_abs * 111.32*1000;
        theta_calc = theta_recalc(rawData);
        % offsets (no angle offset):
        dy = 0.85;
        dx = 0.34;
        % transformation
        X_abs = X_abs + cos(theta_calc) * dx - sin(theta_calc) * dy;
        Y_abs = Y_abs + sin(theta_calc) * dx + cos(theta_calc) * dy;
        % retransform to coordinate angles:
        rawData.LatPos_abs = Y_abs/(111.32*1000);
        rawData.LongPos_abs = X_abs ./ (40075000 * cos(rawData.LatPos_abs*pi()/180) / 360);
        
        save (fullfile(config.root,strcat(matFiles(fileID).name,"_corrected.mat")),'rawData');
    end

    clear data data_table;
    
end