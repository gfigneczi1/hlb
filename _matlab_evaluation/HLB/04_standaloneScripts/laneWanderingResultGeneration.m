close all;
clear;

resultPath = "C:\git\KDP_Igneczi\publikációk\LaneWandering\data\mpcOptimization\resultsSuccessful";

matFiles = dir(fullfile(resultPath, "*.mat"));

for i=1:length(matFiles)
    drID = str2num(matFiles(i).name(7:9));
    load(fullfile(matFiles(i).folder, matFiles(i).name));

    for j=1:length(data)
        if(data(j).e == 'nan')
            e(j) = nan;
        else
            e(j) = data(j).e;
        end
    end
    e_drivers{i} = e(~isnan(e));
    clear e;
end

for i=1:length(e_drivers)
    s_drivers(1:2,i) = [mean(e_drivers{i}); std(e_drivers{i})];
end
disp(s_drivers);