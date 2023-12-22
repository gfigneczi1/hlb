clear;
input_folder =  'C:\git\kdp_hlb_evalframework\_temp';
output_folder = "C:\git\kdp_hlb_evalframework\_temp";

matFiles = dir(fullfile(input_folder,'*.mat'));

for fileID = 1:length(matFiles)
    rawData = load(fullfile(matFiles(fileID).folder,matFiles(fileID).name));
    if (isfield(rawData,'segment'))
        rawData = rawData.segment;
    elseif (isfield(rawData, 'rawData'))
        rawData = rawData.rawData;
    end
    % first adding orientation and curvature based on single track
    % 7500 sample movmean (dt = 10ms, 7500 = 75sec moving window)
    % 0.7 confidence threshold
    x = [0.000, 10.000, 50.000, 85.000, 100.000];
    y = [0.001, 0.100, 0.400, 1.000, 1.000];
    rawData.c1_dx = interp1(x,y, rawData.c1_dx);
    
    doubleTrackAvailable = movmean(rawData.cmass_02, 7500)>0.75; % based on right lane edge
    
    curve_norm_const = 0.00001525878906;
    c1 = ((rawData.c1_orientation .* rawData.c1_dx) + rawData.c2_orientation.*rawData.c2_dx)./ (rawData.c1_dx + rawData.c2_dx);
    c2 = (curve_norm_const * rawData.c1_curvature .* rawData.c1_dx + curve_norm_const * rawData.c2_curvature .* rawData.c2_dx) ./ (rawData.c1_dx + rawData.c2_dx);
    rawData.c1 = zeros(size(rawData.c1_dx));
    rawData.c2 = zeros(size(rawData.c1_dx));
    rawData.c1(doubleTrackAvailable==1) = c1(doubleTrackAvailable==1);
    rawData.c2(doubleTrackAvailable==1) = c2(doubleTrackAvailable==1);
    clear c1 c2
    c1 = (rawData.c1_orientation .* rawData.c1_dx + rawData.c1_orientation.*rawData.c1_dx)./ (rawData.c1_dx + rawData.c1_dx);
    c2 = (curve_norm_const * rawData.c1_curvature .* rawData.c1_dx + curve_norm_const * rawData.c1_curvature .* rawData.c1_dx) ./ (rawData.c1_dx + rawData.c1_dx);
    rawData.c1(doubleTrackAvailable==0) = c1(doubleTrackAvailable==0);
    rawData.c2(doubleTrackAvailable==0) = c2(doubleTrackAvailable==0);
    % weak filtering for connecting issues
    rawData.c1 = movmean(rawData.c1,5);
    rawData.c2 = movmean(rawData.c2,5);
    
    % adding lane distances (c0 values for left and right lane)
    % calculating lane width with strong filtering, wherever the double
    % tracks are available
    estimatedTrackWidth = zeros(size(rawData.c1_dx));
    estimatedTrackWidth(doubleTrackAvailable==1) = movmean(rawData.c01_left(doubleTrackAvailable==1) - rawData.c01_right(doubleTrackAvailable==1), 100000);
    % now filling in the missing parts
    latestTrackWidth = 0; % initialize at zero
    for i=1:length(estimatedTrackWidth)
        if (estimatedTrackWidth(i) > 0)
            latestTrackWidth = estimatedTrackWidth(i);
        else
            estimatedTrackWidth(i) = latestTrackWidth;
        end
    end
    % reconstructing c01_right whenever it is not available
    rawData.c01_right(doubleTrackAvailable==0) = rawData.c01_left(doubleTrackAvailable==0) - estimatedTrackWidth(doubleTrackAvailable==0);
    
    rawData.yawRateESP_sign(rawData.yawRateESP_sign<0) = 1;
    rawData.yawRateESP_sign(rawData.yawRateESP_sign==0) = -1;
    rawData.SteeringAngle = rawData.SteeringAngle * pi() / 180;
    rawData.yawRateESP = rawData.yawRateESP_deg_sec * pi() / 180 .* rawData.yawRateESP_sign;
    rawData = rmfield(rawData,'yawRateESP_sign');
    rawData = rmfield(rawData,'yawRateESP_deg_sec');
    rawData = rmfield(rawData,'c1_dx');
    rawData = rmfield(rawData,'c2_dx');
    rawData = rmfield(rawData,'c1_orientation');
    rawData = rmfield(rawData,'c2_orientation');
    rawData = rmfield(rawData,'c1_curvature');
    rawData = rmfield(rawData,'c2_curvature');
    save(fullfile(output_folder,strcat(matFiles(fileID).name(1:end-4),'_doubleTrack.mat')),'rawData');
end