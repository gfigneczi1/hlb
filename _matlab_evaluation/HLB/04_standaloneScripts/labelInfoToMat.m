clear;
close all;

input_folder = "../../_temp";

% load mat files
matFiles = dir(fullfile(input_folder,'/*.mat'));
mergedData = struct();
for fileID = 1:length(matFiles)
    rawData = load(fullfile(matFiles(fileID).folder,matFiles(fileID).name));
    if (isfield(rawData, "segment"))
        rawData = rawData.segment;
    end
    if (isfield(rawData, "rawData"))
        rawData = rawData.rawData;
    end
    drID = str2num(matFiles(fileID).name(3:5));
    rawData.drID = ones(1, size(rawData.q_T0,2)) * drID;
    fn = fieldnames(rawData);
    if (isempty(fieldnames(mergedData)))
        for i=1:length(fn)
            mergedData.(fn{i}) = rawData.(fn{i});
        end
    else
        for i=1:length(fn)
            if (isfield(mergedData, fn{i}))
                 mergedData.(fn{i}) = [ mergedData.(fn{i}) rawData.(fn{i})];
            end
        end
    end
    
end

save(fullfile(matFiles(fileID).folder,'mergedData.mat'),'-struct', "mergedData");
