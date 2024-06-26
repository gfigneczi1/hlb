function outputStruct = functional_nodePointDefinitionElementary(segment,config, options)
% This function is responsible for finding the optimum node points for one
% given measurement (this is what 'elementary' refers to).
% 4 parts are included:
% - 1: sweep simulation over the parameter space and checking the cost
% surface
% - 2: finding right init value of the node point distances
% - 3: optimizing the node point distances with high steps through the
% entire route and yielding an array of optimum parameters.
% - 4: resimulating the measurement with the optimum parameters.

% @copyright (C) 2022 Robert Bosch GmbH.
% The reproduction, distribution and utilization of this file as
% well as the communication of its contents to others without express
% authorization is prohibited. Offenders will be held liable for the
% payment of damages. All rights reserved in the event of the grant
% of a patent, utility model or design.
% @version 1.0

global startIdx endIdx P3in lookAheadDistance opt segment_in config_in segment_m indexes
% due to MATLAB warning the inputs cannot be global
segment_in = segment;
config_in = config;

%% Options and generic variables
% OBSOLETE: ONLY NEEDED WHEN THREE NODE POINTS ARE USED AND THE SWEEP
% SIMULATION IS ENABLED!!
P1_vector = options.parameters.P1_vector; 
P2_vector = options.parameters.P2_vector;
P3_vector = options.parameters.P3_vector;
startIdx = options.startIdxSurface;
endIdx = options.endIdxSurface;
M = options.independentInitialValues;
stepSize = options.stepSizeLongRange;
P0 = options.parameters.P0;
lookAheadDistance = 250;
opt = 0;

N = round((length(segment_in.X_abs)- 500 - startIdx)/stepSize);
dP = (200-10)/(M-1);

[segment_in, segment_m, indexes] = prepareInputForPlanner(segment);

%% Parameter dependency check
% IMPORTANT: only use with 3 or less node points!
disp("Evaluator - parameter dependency check");
if (options.parameterDependency == "true")
    for i=1:length(P1_vector)
        for j=1:length(P2_vector)
            P30 = [P1_vector(i); P2_vector(j); max(10,250 - P1_vector(i)-P2_vector(j))];
            P30 = P30/250; % Normalization
            fval = functional_optimizeNodePointDistances(P30);
            f_sweep_array_sum_3(i,j) = min(3,fval);
        end
    end
    for i=1:length(P2_vector)
        for j=1:length(P3_vector)
            P30 = [max(10,250 - P1_vector(i)-P2_vector(j)); P1_vector(i); P2_vector(j)];
            P30 = P30/250; % Normalization
            fval = functional_optimizeNodePointDistances(P30);
            f_sweep_array_sum_1(i,j) = min(3,fval);
        end
    end
    for i=1:length(P1_vector)
        for j=1:length(P3_vector)
            P30 = [P1_vector(i); P2_vector(j); max(10,250 - P1_vector(i)-P2_vector(j))];
            P30 = P30/250; % Normalization
            fval = functional_optimizeNodePointDistances(P30);
            f_sweep_array_sum_2(i,j) = min(3,fval);
        end
    end
else
    f_sweep_array_sum_1 = zeros(50,50);
    f_sweep_array_sum_2 = zeros(50,50);
    f_sweep_array_sum_3 = zeros(50,50);
end

% init values
P30 = linspace(options.parameters.mindP, options.parameters.maxDist, options.parameters.Np);
P30 = [options.parameters.mindP, diff(P30)];    
P30 = P30/options.parameters.maxDist; % Normalization

%% Fitting path on N node points
fval = -1;
for j=1:options.parameters.Np
    if (j==1)
        P3in = 1; % normalized with maxDist
    else
        P3in = linspace(options.parameters.mindP, options.parameters.maxDist, j) / options.parameters.maxDist;
        P3in = [options.parameters.mindP, diff(P3in)]; 
    end
    startIdx = options.startIdxInitvalues;
    endIdx = startIdx;
    for i=1:N
        [traj, ref, cor, ~, ~, ~, ~, GT_U, GT_Y]  = functional_driverModelWithPlanner(zeros(j*j,1), 0);
        if (isempty(traj))
        else
            trajResamp = spline(traj(:,1), traj(:,2), cor(:,1));
            refResamp = spline(ref(:,1), ref(:,2), cor(:,1));
            theta = atan(diff(cor(:,2))./diff(cor(:,1))); theta = [theta; theta(end)]; theta = movmean(theta,20);
            refOffset = (refResamp-cor(:,2)).*cos(theta);
            trajOffset = (trajResamp-cor(:,2)).*cos(theta);
            fval = evaluation_lossCalculationTrajectoryPlanner(trajOffset,refOffset,GT_U, GT_Y, -1, P3in);
        end
        f_array_(j,i) = fval;
        startIdx = 1000+i*stepSize;
        endIdx = startIdx; 
    end
end

outputStruct.nodePointCheckFvals = f_array_;

if (options.runGlobalOptimization)
    P3in = [mean(x_array(:,1)) mean(x_array(:,2)) mean(x_array(:,3))];
    
    %% Global optimization - with the right init values
    startIdx = options.startIdxInitvalues;
    endIdx = startIdx;
    for i=1:N
        [x,fval] = fminsearch(@functional_optimizeNodePointDistances,P30);
        f_array_(i) = fval;
        x_array(i,:) = x;
        c2_array(i) = mean(segment_in.curvatureVector(startIdx:startIdx+250));
        d_array(i,:) = (segment_in.X_abs(startIdx:startIdx+250).^2+segment_in.Y_abs(startIdx:startIdx+250).^2).^0.5;
        startIdx = 1000+i*stepSize;
        endIdx = startIdx; 
    end
    P3in = [mean(x_array(:,1)) mean(x_array(:,2)) mean(x_array(:,3))];
    
    %% Resimulation of the measurement with the optimum parameters
    disp("Evaluator - evaluating with optimal parameters");
    [traj, ref, cor, corLeft, corRight, GT_U, GT_Y] = functional_driverModelWithPlanner(P0, 0);
    f = evaluation_lossCalculationTrajectoryPlanner(traj(:,:),ref(:,:),GT_U, GT_Y, opt, P3in);
    disp(strcat('loss with resimulated global parameters:'));
    disp(f);
    
    %% Saving the output signals
    outputStruct.f = f;
    outputStruct.optimization.x = x_array;
    outputStruct.optimization.f = f_array_;
    outputStruct.optimization.kappa = c2_array;
    outputStruct.optimization.d = d_array;
    outputStruct.traj = traj;
    outputStruct.ref = ref;
    outputStruct.cor = cor;
    outputStruct.corLeft = corLeft;
    outputStruct.corRight = corRight;
    outputStruct.P = P0;
    outputStruct.P3in = P3in;
    outputStruct.initvalueCosts = f_big_array;
    outputStruct.paramDependency.surf_1 = f_sweep_array_sum_1;
    outputStruct.paramDependency.surf_2 = f_sweep_array_sum_2;
    outputStruct.paramDependency.surf_3 = f_sweep_array_sum_3;
    outputStruct.paramDependency.P1 = P1_vector;
    outputStruct.paramDependency.P2 = P2_vector;
    outputStruct.paramDependency.P3 = P3_vector;
    outputStruct.anchors.U = GT_U;
    outputStruct.anchors.Y = GT_Y;
end
end

