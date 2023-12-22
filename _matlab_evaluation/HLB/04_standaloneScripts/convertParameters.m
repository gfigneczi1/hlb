function convertParameters()
    paramFiles = dir(fullfile('C:\git\kdp_hlb_evalframework\_temp','parameters*.mat'));
    for fileID = 1:length(paramFiles)
            parameterRawData = load(fullfile(paramFiles(fileID).folder,paramFiles(fileID).name));
    end
end

