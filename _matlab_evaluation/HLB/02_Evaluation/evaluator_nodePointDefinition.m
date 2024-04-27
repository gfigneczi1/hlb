function evaluator_nodePointDefinition(segments,config)
% This evaluator is responsible to produce the optimum node point distance
% set for the given driver(s). The optimization happens through various
% costs defined in evaluation_lossCalculationTrajectoryPlanner.
% One ore more measurements can be used.

% Segmentor: segmentor_driverModel.

% @copyright (C) 2022 Robert Bosch GmbH.
% The reproduction, distribution and utilization of this file as
% well as the communication of its contents to others without express
% authorization is prohibited. Offenders will be held liable for the
% payment of damages. All rights reserved in the event of the grant
% of a patent, utility model or design.
% @version 1.0

global P3in segment_m startIdx endIdx model_ver options

model_ver = 0; % ground truth

options.parameters.P1_vector = linspace(0,240,50);
options.parameters.P2_vector = linspace(0,240,50);
options.parameters.P3_vector = linspace(0,240,50);
options.parameters.Np = 15;
options.parameters.mindP = 10;
options.parameters.maxDist = 150;
options.startIdxSurface = 3000;
options.endIdxSurface = 3000;
options.startIdxInitvalues = 1000;
options.startIdxGlobal = 1000;
options.independentInitialValues = 21;
options.stepSizeLongRange = 200;
options.parameters.P0 = zeros(21,1); % LDM parameter set
options.parameters.P3in = [0; 28.576; 98.268]; % node point parameter set
options.parameterDependency = "false";
options.runGlobalOptimization = false;

%% Node point optimization
% Calling the node point definition for different segments
% Config and options are the same
for i=1:length(segments.segments)
    disp(strcat("Evaluating file",{' '},num2str(i),'/',num2str(length(segments.segments))));
    segment = segments.segments(i).segment;
    outputStruct{i} = functional_nodePointDefinitionElementary(segment,config, options);
    if (options.runGlobalOptimization)
        p1(i) = max(10/250,outputStruct{i}.P3in(1));
        p2(i) = outputStruct{i}.P3in(2);
        p3(i) = outputStruct{i}.P3in(3);
        disp(strcat('Node point distances for', {' '}, segments.segments(i).name, {' '},'are:'));
        disp([(p1(i)*250) (p1(i)+p2(i))*250 (p1(i)+p2(i)+p3(i))*250]);
    end
end

plot_nodePointNumberAnalysis(outputStruct, options, config);

if (options.runGlobalOptimization)
%% Global optimum parameters - together for all drivers
% Calculate the characteristic distance set and adding the results for all
% drivers
P3Global = [mean(p1) mean(p2) mean(p3)];
disp('Global node point distances are:');
disp([P3Global(1) P3Global(1)+P3Global(2) P3Global(1)+P3Global(2)+P3Global(3)]*250);

for i=1:length(segments.segments)
    segment = segments.segments(i).segment;
    [~, segment_m, ~] = prepareInputForPlanner(segment);
    
    P0 = zeros(21,1);
    P3in = P3Global;
    
    startIdx = options.startIdxGlobal;
    N = round((length(segment.X_abs)- 500 - startIdx)/options.stepSizeLongRange);
    endIdx = 1000+N*options.stepSizeLongRange; 
    
    [traj, ref, ~, ~, ~, GT_U, GT_Y] = functional_driverModelWithPlanner(P0, 0);
    f = evaluation_lossCalculationTrajectoryPlanner(traj(:,:),ref(:,:),GT_U, GT_Y, 2, P3in);

    outputStruct{i}.global.P3Global = P3Global;
    outputStruct{i}.global.f = f;
    outputStruct{i}.global.traj = traj;
    outputStruct{i}.global.ref = ref;
    outputStruct{i}.name = segments.segments(i).name;
end
end

%% Evaluation
% saving raw evaluation data
save(fullfile(config.root,'plots','outputRawStruct.mat'),'outputStruct');

%% Plotting
set(0,'DefaultFigureVisible','off');
plot_nodePointDefinition(outputStruct,  options, config);
set(0,'DefaultFigureVisible','on');

end


