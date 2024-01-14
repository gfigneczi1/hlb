close all; clear; 

%path = "C:\git\KDP_Igneczi_new\KDP_Igneczi\publikációk\LaneWandering\data\offsetSeparation";
path = "C:\git\hlb\_temp\plots\frequencies";

matFiles = dir(fullfile(path, "*.mat"));

frequencies = [linspace(0.016, 0.098, 5), linspace(0.114, 1, 10)];

for i=1:length(matFiles)
    matFile = fullfile(matFiles(i).folder, matFiles(i).name);
    data = load(matFile);
    for j=1:length(data.data)
        coverage(i,j) = data.data(j).coverage;
        sigma_c(i,j) = data.data(j).lengths(2)';
    end
end

w_c = 0.5;
w_s = 0.5;
cost = w_c* (1-coverage) + w_s*sigma_c/6;