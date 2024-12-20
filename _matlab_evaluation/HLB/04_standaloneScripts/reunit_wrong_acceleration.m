clear;

matFiles = dir(fullfile("C:\git\KDP\HLB_for_AV_MotionControl\02_Results\MassMeasurements\Dr006\20241120", "*withTraffic.mat"));

for i=1:length(matFiles)
    data = load(fullfile(matFiles(i).folder, matFiles(i).name));
    if (isfield(data, 'segment'))
        data = data.segment;
    end
    data.AccelerationX = data.AccelerationX*9.81;
    data.Acceleration_Y = data.Acceleration_Y*9.81;
    save(fullfile(matFiles(i).folder,matFiles(i).name), '-struct', 'data');
end
    