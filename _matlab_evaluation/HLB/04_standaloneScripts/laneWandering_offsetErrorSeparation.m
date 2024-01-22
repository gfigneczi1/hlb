close all; clear; 

%path = "C:\git\KDP_Igneczi_new\KDP_Igneczi\publikációk\LaneWandering\data\offsetSeparation";
path = "C:\database\LaneWand\frequencies";

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

f = figure();

set(0,'defaultAxesFontSize',14);
set(f,'defaulttextInterpreter','latex') ;
set(f, 'defaultAxesTickLabelInterpreter','latex');  
set(f, 'defaultLegendInterpreter','latex');

minCov = min(coverage);
maxCov = max(coverage);
meanCov = mean(coverage);

minSig = min(sigma_c);
maxSig = max(sigma_c);
meanSig = mean(sigma_c);

fill([frequencies flip(frequencies)], [1-maxCov, flip(1-minCov)], 'y', 'FaceAlpha',0.5, 'DisplayName', '$1-coverage$ range'); hold on;
plot(frequencies, 1-meanCov, 'k', 'LineWidth', 2, 'DisplayName', '$1-coverage$');

fill([frequencies flip(frequencies)], [minSig/6, flip(maxSig/6)], 'b', 'FaceAlpha',0.5, 'DisplayName', '$\sigma_c$ range'); hold on;
plot(frequencies, meanSig/6, 'k', 'LineWidth', 2, 'DisplayName', '$\frac{\sigma_c}{\sigma_c^{max}}$'); 

legend;

grid on;
