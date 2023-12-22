close all;
clear all;

path = "C:\git\KDP_Igneczi\publikációk\LaneWandering\data\mpcOptimization\resimulation_MWM_MOE";
qualityMap = "C:\git\KDP_Igneczi\publikációk\LaneWandering\data\mpcOptimization\quality.xlsx";

q = table2array(readtable(qualityMap));
q = q(2:end,2:end);

matFiles = dir(fullfile(path, "\*.mat"));

f = figure();

for i=1:length(matFiles)
    data = load(fullfile(matFiles(i).folder, matFiles(i).name));
    qualityVector = q(i,:);
    e = [];
    for j=1:size(data.data,2)
        if (qualityVector(j)==1)
            e = [e; data.data(j).e];
        end
    end
    subplot(2,1,1);
    plot(e, "DisplayName",num2str(i), 'Marker', 'x');
    hold on;
    grid on;
    ratio(i) = numel(find(e<=0.04))/length(e);
end

title("RMS of the Euler distance between planned and original path");

xlabel('snippet ID'); ylabel('RMS of ED(m)');
ylim([0,0.5]);

subplot(2,1,2);
plot(ratio, 'bo');
grid on;
xlabel('Driver ID');
ylabel('ratio(-)');
ylim([0,1]);
xticks([1:1:length(ratio)]);
