close all; clear; 

path = "C:\git\KDP_Igneczi_new\KDP_Igneczi\publikációk\LaneWandering\data\offsetSeparation";

matFiles = dir(fullfile(path, "*.mat"));

for i=1:length(matFiles)
    matFile = fullfile(matFiles(i).folder, matFiles(i).name);
    data = load(matFile);
    for j=1:length(data.data)
        coverage(i,j) = data.data(j).coverage;
        sigma_c(i,j) = data.data(j).lengths(2)';
    end
end