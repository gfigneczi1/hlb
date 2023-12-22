function evaluator_driverModelSimulation(segmentInAndOut, config)
global segment_m lookAheadDistance startIdx endIdx GT_U GT_theta GT_Y indexes anchorPoints_array Poptim
global P3in;
P3in = [0; 26.59; 96.91];
%P3in = [89.0071; 117.6216;   51.7895];

model_input_struct = segmentInAndOut.input;
model_output = segmentInAndOut.output;

lookAheadDistance = 150;

segment = struct2table(model_input_struct);
clear segmentStruct;

for j=1:size(segment.c2,1)
    corridor(j,1:2) = pos_tf2GPS(segment.X_abs(j),segment.Y_abs(j),segment.theta_calc(j),0.5*(segment.c01_left(j)+segment.c01_right(j)));
end

segment.corrX = corridor(:,1);
segment.corrY = corridor(:,2);
segment.orientationVector = atan(pt1_filter(segment.c1,1)) + segment.theta_calc;
segment.curvatureVector = movmean(segment.c2,100);

% SAVING out field indices
varnames = segment.Properties.VariableNames;
[~, indexes.X_abs] = ismember('X_abs', varnames);
[~, indexes.Y_abs] = ismember('Y_abs', varnames);
[~, indexes.theta_calc] = ismember('theta_calc', varnames);
[~, indexes.c1] = ismember('c1',varnames);
[~, indexes.c2] = ismember('c2',varnames);
[~, indexes.orientationVector] = ismember('orientationVector',varnames);
[~, indexes.curvatureVector] = ismember('curvatureVector',varnames);
[~, indexes.corrX] = ismember('corrX',varnames);
[~, indexes.corrY] = ismember('corrY',varnames);
[~, indexes.velocityX] = ismember('VelocityX_ESP',varnames);
[~, indexes.q_T0] = ismember('q_T0',varnames);
[~, indexes.yawRate] = ismember('yawRateESP',varnames);

segment_m = table2array(segment);
segment_m(:,end+1) = segment.yawRateESP/0.2222.*movmean(segment.SteeringTorque/3,100);

startIdx = 2;
endIdx = size(segment,1);

P = zeros(21,1);
[f, traj, ref, cor] = objectivefcn1(P, 1, 0, "");
evaluation_ValidateSimulationResults(traj, model_output);
end