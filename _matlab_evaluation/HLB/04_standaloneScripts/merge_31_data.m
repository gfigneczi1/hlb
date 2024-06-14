matFiles = dir("C:\git\KDP\HLB_for_AV_MotionControl\02_Results\MassMeasurements\Dr021\standardized\*r31_withTraffic.mat");

outputName = "Dr021_2023-04-26_r31_withTraffic.mat";
for i=1:length(matFiles)
    data = load(fullfile(matFiles(i).folder,matFiles(i).name));
    
    fn = fieldnames(data);
    if (i==1)
        dataMerged = data;
    else
        for j=1:length(fn)
            dataMerged.(fn{j}) = [dataMerged.(fn{j}) data.(fn{j})];
        end
    end
end

save(fullfile(matFiles(i).folder, outputName), '-struct', 'dataMerged');