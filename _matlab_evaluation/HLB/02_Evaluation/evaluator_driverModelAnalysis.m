function evaluator_driverModelAnalysis(segments,config)

% This evaluator is responsible for generating correlation plots of a human
% driver. The following plots are generated:
% - node point offsets vs. curvatures
% - node point offsets vs. curvature gradients
% - trajectory plots (maps)
% The script is prepared to handle multiple segments.

% Segmentor: segmentor_driverModel.

% @copyright (C) 2022 Robert Bosch GmbH.
% The reproduction, distribution and utilization of this file as
% well as the communication of its contents to others without express
% authorization is prohibited. Offenders will be held liable for the
% payment of damages. All rights reserved in the event of the grant
% of a patent, utility model or design.
% @version 1.0
% @date 2022-12-07

global lookAheadDistance P3in startIdx endIdx segment_m indexes model_ver

% Init definitions
GT_U_array = []; GT_Y_array = [];
P3in = [0; 26.59; 96.91];
%P3in = [61.3546; 114.1398;   53.1781];
P3in = P3in/250;
P0 = zeros(21,1); % this is a placeholder only, as the parameter does not play any role when using driver data directly
lookAheadDistance = 150; % in meters
startIdx = 100;
stepSize = 300;

model_ver = 0;

for fileID=1:size(segments.segments,2)
    segment = segments.segments(fileID).segment;
    name = segments.segments(fileID).name;
    signalStatus = segments.segments(fileID).signalStatus;
    cutInfo = segments.segments(fileID).cutInfo;
    
    try
    
        [segment_in, segment_m, indexes] = prepareInputForPlanner(segment);

        if (length(segment_in.X_abs) < 1000)
            startIdx = 2;
            endIdx = length(segment_in.X_abs)-startIdx;
            N = 1;
        else
            N = round((length(segment_in.X_abs) - startIdx)/stepSize);
            endIdx = min(startIdx+N*stepSize,length(segment_in.X_abs)-startIdx);
        end
        if (N > 0)
            [traj, ref, cor, corLeft, corRight, ~, curv, GT_U, GT_Y, validPoints, replan_array, GT_U_array, GT_Y_array, trajFull, refFull, corFull] = functional_driverModelWithPlanner(P0, 0);

            % measurement based evaluation
            if (~isempty(traj))
                plot_trajectoryComparison(config,ref,traj,cor,corLeft,corRight,strcat('undisturbed_',name));
                plot_correlationPlotsLDM(config,GT_U,GT_Y,strcat('undisturbed_',name));
            end
            
            %[traj, ref, cor, corLeft, corRight, ~, ~, GT_U, GT_Y] = functional_driverModelWithPlanner(P0, 3);

            % measurement based evaluation
            if (~isempty(traj))
                plot_trajectoryComparison(config,ref,traj,cor,corLeft,corRight,strcat('disturbed_',name));
                plot_correlationPlotsLDM(config,GT_U,GT_Y,strcat('disturbed_',name));
            end
        end
        plot_cutInfo(config, segments.segments(fileID).uncut, cutInfo, name);
        
        plot_driverModelTrajectoryHistogram([], traj(abs(curv)>2.5e-4,:),cor(abs(curv)>2.5e-4,:), strcat(name,"_deviationHistogramHumanInCurve"));
        plot_driverModelTrajectoryHistogram([], traj(abs(curv)<=2.5e-4,:),cor(abs(curv)<=2.5e-4,:), strcat(name,"_deviationHistogramHumanInStraight"));

        parameters_out(fileID).U = GT_U;
        parameters_out(fileID).dY = GT_Y(2:end,:);
        parameters_out(fileID).drID = str2num(name(3:5));

        clear segment;
    catch e
        disp("The following file could not be evaluated:");
        disp(segments.segments(fileID).name);
        disp(e);
    end
end

save(fullfile(config.root,"plots", 'parametersAllDrivers.mat'), 'parameters_out');

% all measurements evaluation
%plot_correlationPlotsLDM(config,GT_U_array,GT_Y_array,"all");
    
end


