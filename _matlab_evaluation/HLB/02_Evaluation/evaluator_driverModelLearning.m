function evaluator_driverModelLearning(segments,config)
% This evaluator is responsible for generating doing the learning 
% parameters of LDM based on human driving data. The following plots are
% generated:
% - curvature deviation in curves
% - position deviation in curves
% - trajectory shapes in curves
% - histrograms on position deviation
% A report is generated based on the KPI-s defined for the model.
% The script is prepared to handle multiple segments. Learning is done
% separately for the files.

% Segmentor: segmentor_driverModel.

% @copyright (C) 2022 Robert Bosch GmbH.
% The reproduction, distribution and utilization of this file as
% well as the communication of its contents to others without express
% authorization is prohibited. Offenders will be held liable for the
% payment of damages. All rights reserved in the event of the grant
% of a patent, utility model or design.
% @version 1.0
% @date 2022-12-08

global lookAheadDistance P3in segment_m indexes startIdx endIdx model_ver GT_U_GP alpha Kxx GT_Y_GP

lookAheadDistance = 150; % in meters
P0 = zeros(21,1); % this is a placeholder only, as the parameter does not play any role when using driver data directly
startIdx = 2;
stepSize = 300;
%P3in = [0; 26.59; 96.91]; % node point distances from previous optimization
P3in = [20; 40; 80];
P3in = P3in/250;
model_ver = 0;

for fileID=1:size(segments.segments,2)
    segment = segments.segments(fileID).segment;
    name = segments.segments(fileID).name;
   
    try
        [segment_in, segment_m, indexes] = prepareInputForPlanner(segment);

        if (length(segment_in.X_abs) < 1000)
            startIdx = 2;
            endIdx = length(segment.X_abs)-startIdx;
            N = 1;
        else
            N = round((length(segment_in.X_abs) - startIdx)/stepSize);
            endIdx = min(startIdx+N*stepSize,length(segment_in.X_abs)-startIdx);
        end

        %% GROUND TRUTH - Learning basis
        [traj, ref, cor, ~, ~, orient, curv, GT_U, GT_Y, validPoints, ~, GT_U_array, GT_Y_array, ~, ~, ~] = functional_driverModelWithPlanner(P0, 0);
        driverModelOutput.traj = traj;
        driverModelOutput.ref = ref;
        driverModelOutput.cor = cor;
        driverModelOutput.orient = orient;
        driverModelOutput.curv = curv;
        driverModelOutput.cor = cor;
        driverModelOutput.GT_U = GT_U_array;
        driverModelOutput.GT_Y = GT_Y_array;
        GT_U_GP = GT_U;
        GT_Y_GP = GT_Y;
        
        
        % Driver evaluation (ground truth)
        KPI = evaluate_KPIs_LDM(segment, driverModelOutput, strcat(name,'_Human_'), config,[]);
        plot_correlationPlotsLDM(config,GT_U,GT_Y,strcat(name,'_Human_'));
        
        GT_U_driver = GT_U;
        GT_Y_driver = GT_Y;
        
        
        % Plotting
        % Plot generic data
        plot_driverModelLearning(driverModelOutput,KPI, config, segment, strcat(name,'_Human'), []);
        
        %% CURVE based model verification
        model_ver = 1;
        if (model_ver <=1)
            % learning the Parameters of given driver
            if (model_ver == 0)
                Poptim = functional_driverModelLearning(GT_U_driver,GT_Y_driver,2);
            elseif (model_ver == 1)
                Poptim = functional_driverModelLearning(GT_U_driver,GT_Y_driver,7);
            end
            parameters.(name(1:5)) = Poptim;
        elseif (model_ver == 2)
            global GP_mode;
            Kxx = [];
            GP_mode = 0; % mean behavior
            %% Hyper parameters and kernel selection
            sigma_es=0.15;
            sigma_ker=0.2;
            Poptim = zeros(21,1);
            kernel = @(x1,x2)  exp(-(x1-x2)'*(x1-x2)/(2*sigma_ker));
            N = size(GT_U_GP,2);

            for k=1:N
                for j=1:N
                    Kxx(k,j)=kernel(GT_U_GP(1:3,k),GT_U_GP(1:3,j));
                end
            end

            alpha=(Kxx+sigma_es*eye(N))\GT_Y_GP(2:end,:)';
        elseif (model_ver == 3)
            %% Define GP 
            meanfunc = [];          % Start with a zero mean prior
            covfunc = @covSEiso;    % Squared Exponental covariance function
            Poptim = zeros(21,1);
            global hyp_opt1 hyp_opt2 hyp_opt3 

            %% ID problem
            likfunc = @likGauss;    % Gaussian likelihood
            hyp = struct('mean', [], 'cov', [0 0], 'lik', -1); % initalize hyper-parameter structure associated with the mean, covariance an likelihood functions
            hyp_opt1 = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, GT_U_GP(1:3,:), GT_Y_GP(2,:)); % Optimize the marginal likelihood
            hyp_opt2 = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, GT_U_GP(1:3,:), GT_Y_GP(3,:)); % Optimize the marginal likelihood
            hyp_opt3 = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, GT_U_GP(1:3,:), GT_Y_GP(4,:)); % Optimize the marginal likelihood
        end
        %% Model validation
        clear GT_U GT_Y KPI
        if (model_ver <=1)
            [traj, ref, cor, ~, ~, orient, curv, GT_U, GT_Y,  ~, ~] = functional_driverModelWithPlanner(Poptim, 1);
            driverModelOutput.traj = traj;
            driverModelOutput.ref = ref;
            driverModelOutput.cor = cor;
            driverModelOutput.orient = orient;
            driverModelOutput.curv = curv;
            driverModelOutput.cor = cor;
            try
                % curve based evaluation
                KPI{fileID} = evaluate_KPIs_LDM(segment, driverModelOutput, strcat(name,'_Model_'), config, []);
                KPI{fileID}.name = name;

                %% Plotting
                % Plot generic data
                plot_driverModelLearning(driverModelOutput,KPI{fileID}, config, segment, strcat(name,'_Model_'), []);
                % Plot correlation data
                %plot_correlationPlotsLDM(config,GT_U,GT_Y,strcat(name,'_Model_', num2str(i)));

                plot_driverCorrelationWithModel(config,GT_U_driver, GT_Y_driver, GT_U, GT_Y, strcat(name,'_Comparison_', []));
            catch e
                disp('Plotting and KPI evaluation is unsucessful for the following file in the following curve:');
                disp(name);
            end
        elseif (model_ver >= 2)
            [traj, ref, cor, ~, ~, orient, curv, GT_U, GT_Y,  ~, ~] = functional_driverModelWithPlanner(Poptim, 1);
            driverModelOutput.traj = traj;
            driverModelOutput.ref = ref;
            driverModelOutput.cor = cor;
            driverModelOutput.orient = orient;
            driverModelOutput.curv = curv;
            driverModelOutput.cor = cor;
            [curves] = cutCurves(driverModelOutput);
            for i=1:length(curves)
                try
                    % curve based evaluation
                    KPI{fileID} = evaluate_KPIs_LDM(segment, driverModelOutput, strcat(name,'_Model_'), config, i);
                    KPI{fileID}.name = name;

                    %% Plotting
                    % Plot generic data
                    plot_driverModelLearning(driverModelOutput,KPI{fileID}, config, segment, strcat(name,'_Model_'), i);
                    % Plot correlation data
                    %plot_correlationPlotsLDM(config,GT_U,GT_Y,strcat(name,'_Model_', num2str(i)));

                    plot_driverCorrelationWithModel(config,GT_U_driver, GT_Y_driver, GT_U, GT_Y, strcat(name,'_Comparison_', num2str(i)));
                catch
                    disp('Plotting and KPI evaluation is unsucessful for the following file in the following curve:');
                    disp(name);
                    disp(strcat("Curve", num2str(i)));
                end
            end
        end
%         if (model_ver == 2)
%             GP_mode = 1; % upper sigma
%             [trajUpper, ~, ~, ~, ~, ~, ~, GT_U_upper, GT_Y_upper,  ~, ~] = functional_driverModelWithPlanner(Poptim, 1);
%             GP_mode = -1; % upper sigma
%             [trajlower, ~, ~, ~, ~, ~, ~, GT_U_lower, GT_Y_lower,  ~, ~] = functional_driverModelWithPlanner(Poptim, 1);
%         end

    catch e
        disp('Evaluation is unsuccessfuly for the following file:');
        disp(name);
        disp(e);
    end
end

save(fullfile(config.root, 'plots', 'parameters.mat'), 'parameters');