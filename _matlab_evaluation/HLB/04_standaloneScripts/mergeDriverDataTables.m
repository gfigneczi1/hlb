function mergeDriverDataTables(token)
    folder = "C:\git\kdp_hlb_evalframework\_temp\plots";

    files = dir(fullfile(folder,token));

    for i=1:length(files)
        t = xlsread(fullfile(files(i).folder, files(i).name));
        T(i,:) = t;
    end

    xlswrite(fullfile(folder,strcat(files(1).name(1:end-7),'.xlsx')),T);
end