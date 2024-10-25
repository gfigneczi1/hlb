clear;

matFiles = dir(fullfile("C:\git\KDP\HLB_for_AV_MotionControl\02_Results\MassMeasurements\Dr023", "*.mat"));

for i=1:length(matFiles)
    data = load(fullfile(matFiles(i).folder, matFiles(i).name));
    data.AccelerationX_ESP = data.AccelerationX_ESP*9.81;
    data.AccelerationY_ESP = data.AccelerationY_ESP*9.81;
    data.yawRateESP = data.yawRateESP*pi()/180;
    data.yawAngle = data.yawAngle*pi()/180;
    save(fullfile(matFiles(i).folder,matFiles(i).name), '-struct', 'data');
end
    