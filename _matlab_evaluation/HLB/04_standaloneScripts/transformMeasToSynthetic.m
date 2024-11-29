close all; clear; clc;

load('C:\database\KDP_HLB_GP\Dr008_Dr023_input_GP_fixedAcceleration.mat');

for i=1:length(segments.segments)
    segment = segments.segments(i).segment;
    LaneCurvatureGradient = movmean(diff(tan(segment.LaneCurvature))./diff(segment.X_abs),200); LaneCurvatureGradient = [LaneCurvatureGradient; LaneCurvatureGradient(end)];
    input = [segment.OncomingTrafficType ...
        segment.FrontTrafficType ...
        segment.VelocityX ...
        segment.AccelerationX ...
        segment.YawRate ...
        segment.LaneCurvature ...
        LaneCurvatureGradient];

    data.input = input;
    data.egoPosition = [segment.X_abs segment.Y_abs];
    data.laneOrientation = segment.LaneOrientation;
    data.relativeTime = segment.Relative_time;

    save(strcat("syntheticFromRealMeas_", num2str(i), ".mat"), "data");
end