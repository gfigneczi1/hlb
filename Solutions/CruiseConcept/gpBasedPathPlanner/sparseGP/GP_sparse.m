% This function is the cover function to call the GP sub-function for a
% given segment_m. segment_ms is an array of struct containing the driving data
% of multiple drivers. Config is a metadata struct containing information
% about e.g., the plot folders
close all; clear; clc;
addpath(fullfile("..", "library"));
addpath(fullfile("..", "library", "gpml-matlab-master"));
addpath(fullfile("..", "library", "gpml-matlab-master", "cov"));
addpath(fullfile("..", "library", "gpml-matlab-master", "inf"));
addpath(fullfile("..", "library", "gpml-matlab-master", "lik"));
addpath(fullfile("..", "library", "gpml-matlab-master", "mean"));
addpath(fullfile("..", "library", "gpml-matlab-master", "prior"));
addpath(fullfile("..", "library", "gpml-matlab-master", "util"));
load('Dr008_Dr027_input_GP.mat');
if (~exist("plots", 'dir'))
    mkdir("plots");
else
    rmdir("plots", 's');
    mkdir("plots");
end
if (~exist("kpis", 'dir'))
    mkdir("kpis");
else
    rmdir("kpis", 's');
    mkdir("kpis");
end

%% SCRIPT CONFIGURATION
DRIVER_ID = []; % [1, 5, 6, 10, 11, 13]; % possible values: 1,2,3,4,5...etc or []=all drivers simulated
ITERATE_INDUCTION = true;

%% DATA PIPELINE PARAMETERS
MAXIMUM_INPUT_SIZE = inf; % before snipetting and normalization
KERNEL_TYPE = "{'covPPard',3}";
RATIO_OF_TRAIN_VS_TOTAL = 0.7;
MAXIMUM_PREVIEW_DISTANCE = 150; % from within preview information is extracted, unit: m
OUTPUT_SHIFT = linspace(15,MAXIMUM_PREVIEW_DISTANCE,10); %possible values: scalar number or array e.g., linspace(15,MAXIMUM_PREVIEW_DISTANCE,10);
INDUCED_SIZE = 100; % only applicable if ITERATE_INDUCTION is FALSE
inductionSizePool = [10, 20, 50];%, 100, 300, 500, 1000]; % only applicable if ITERATE_INDUCTION is TRUE

variables = ["t_{ot}", "f_{ot}", "v_x", "a_x", "\omega", "\kappa", "d\kappa"];

if (isempty(DRIVER_ID))
    i0 = 9; i1 = length(segments.segments);
    if (i0~=1)
        fprintf("WARNING: DRIVER ID 0 is not zero, simulating only from driver 3!!!\n");
    end
elseif (~isempty(DRIVER_ID) && length(DRIVER_ID) > 1)
    i0 = 1; i1 = 1;
elseif (~isempty(DRIVER_ID) && DRIVER_ID<=length(segments.segments)-2)
    i0 = DRIVER_ID; i1=i0;
else
    i0 = []; i1 = [];
    disp("driver ID is invalid!")
end
if (~isempty(i0) && ~isempty(i1))
    if (~isempty(DRIVER_ID))
        segmentMerged=segments.segments(DRIVER_ID).segment;
    end
    if (length(DRIVER_ID) > 1)
        for i=1:length(DRIVER_ID)
            % loop through segment_ms and concatenate each to the previous one
            if (i==i0)
                segmentMerged=segments.segments(DRIVER_ID(i)).segment;
            else
                segmentMerged=[segmentMerged; segments.segments(DRIVER_ID(i)).segment];
            end
        end
    else
        for i=i0:i1
            % loop through segment_ms and concatenate each to the previous one
            if (i==i0)
                segmentMerged=segments.segments(i).segment;
            else
                segmentMerged=[segmentMerged; segments.segments(i).segment];
            end
        end
    end
    [~, segmentMerged_m, indexesMerged] = prepareInputForPlanner(segmentMerged);
    
    PARAMS.MAXIMUM_INPUT_SIZE = MAXIMUM_INPUT_SIZE;
    PARAMS.OUTPUT_SHIFT = OUTPUT_SHIFT;
    [~, ~, ~, ~, ~, s_out, ~, s_in] = prepareData(segmentMerged_m, indexesMerged, PARAMS);
    
    set(0,'DefaultFigureVisible','off');
    
    % PARAMETERS REMAPPING
    shiftOnOutputSelection = OUTPUT_SHIFT;  %shift the offset in time (positive means shift forward)
    p = RATIO_OF_TRAIN_VS_TOTAL; %percentage of evaluation data from entire dataset
    
    if (ITERATE_INDUCTION)
        indSize0 = 2; indSize1 = 2; %length(inductionSizePool);
    else
        indSize0 = 1; indSize1 = 1;
    end
    for indSize = indSize0:indSize1
        if (ITERATE_INDUCTION)
            INDUCED_SIZE = inductionSizePool(indSize);
        end
        for driverID = i0:i1
            if (length(DRIVER_ID)>1)
                segment_m = segmentMerged_m;
                indexes = indexesMerged;
                name = "mergedDrivers";
            else
                segment = segments.segments(driverID).segment;
                name = segments.segments(driverID).name;
    
                % transforming the struct to array for lower calculation time
                [~, segment_m, indexes] = prepareInputForPlanner(segment);
            end            
            dT = mean(diff(segment_m(:, indexes.Relative_time)));
        
            for shiftID=1:numel(OUTPUT_SHIFT)
                tic;
                dx = segment_m(:, indexes.VelocityX)*dT;
                shiftOnOutput = [1:1:size(segment_m(:,indexes.Relative_time),1)]'+floor(OUTPUT_SHIFT(shiftID)./dx);
                input = [segment_m(:, indexes.OncomingVehicleTimeToPass), ...
                        segment_m(:, indexes.OncomingTrafficType), ...
                        segment_m(:, indexes.FrontTrafficType), ...
                        segment_m(:, indexes.VelocityX), ...
                        movmean(segment_m(:, indexes.AccelerationX),20), ...
                        movmean(segment_m(:, indexes.YawRate),20), ...
                        movmean(segment_m(:, indexes.LaneCurvature), 20), ...
                        movmean(segment_m(:, indexes.c3), 200)];
                input(:,2) = movmean(input(:,2),50);
                input(:,3) = movmean(input(:,3),50);
                
                input = input(:,2:end);
                inputRaw = input; % saved for further analysis
        
                output = zeros(size(input,1),1);
                delta = -segment_m(:, indexes.c0);
                for shiftIDonOutput=1:size(input,1)
                    if (shiftOnOutput(shiftIDonOutput) > size(input,1))
                        break;
                    else
                        output(shiftIDonOutput,1) = delta(shiftOnOutput(shiftIDonOutput));
                    end
                end   
                outputRaw = output; % saved for further analysis
            
                % SHUFFLE
                N = size(input,1);
                shuffledIndeces = randperm(N);
                input = input(shuffledIndeces,:);
                output = output(shuffledIndeces, :);
        
                % LIMIT DATA IF NEEDED
                % this is done before norm and central
                input = input(1:min(size(input,1),MAXIMUM_INPUT_SIZE),:);
                output = output(1:min(size(input,1),MAXIMUM_INPUT_SIZE), :);
        
                % NORM AND CENTRAL
                % According to the global cin and sin, cout and sout
                [~, cout, ~] = normalize(output);
                [~, cin, ~] = normalize(input);
                output = (output-cout)./s_out(shiftID);
                input = (input-cin)./s_in;
                
                % REMOVING SPARE POINTS I.E., OUTPUT HIGHER THAN 2, = 2 sigma
                % variance
                input = input(abs(output(:,1))<=3,:);
                output = output(abs(output(:,1))<=3,:);

                plotCorrelation (input, output, shiftID, driverID, variables);
        
                % SIMULATION
                M = size(input,1);
                % EVALUATION / VALIDATION DATA SELECTION
                estimationData = 1:1:round(p*M);
                validationData = round(p*M)+1:1:M;
                input_estimation = input(estimationData,:);
                output_estimation = output(estimationData,1);
                input_validation = input(validationData,:);
                output_validation = output(validationData,:);
        
                %% Define GP 
                hyp = generateHyperParamStruct(KERNEL_TYPE, size(input,2));
                eval(strcat('covfunc = ',KERNEL_TYPE));    % Squared Exponental covariance function
                meanfunc = [];       % Start with a zero mean prior
                likfunc = @likGauss;    % Gaussian likelihood
                
                xu = input(1:INDUCED_SIZE,:);
                hypInduced = hyp;
                covInduced = {'apxSparse',covfunc,xu};           % inducing points
                inf = @infGaussLik;
                infv  = @(varargin) inf(varargin{:},struct('s',0.0));           % VFE, opt.s = 0
                hypInduced.xu = xu;  % adding induced inputs
                hyp_opt = minimize(hypInduced,@gp,-100, inf,meanfunc,covInduced,likfunc,input,output); 
        
                input_induced= hyp_opt.xu;
                output_induced = gp(hyp_opt,infv,meanfunc,covInduced,likfunc, input, output, input_induced);  
            
                % checking performance degradation
                [estimationInduced, deviationInduced] = gp(hyp_opt, infv, meanfunc, covInduced, likfunc, input_induced, output_induced, input_estimation);
        
                %% KPI-s
                [RMSinduced, NRMS_W_induced, NRMS_M_induced, RMS_DEV_induced] = KPIcalculation(estimationInduced, deviationInduced, output_estimation);
                fprintf("Induction accuracy:\n");
                fprintf("RMS value is: %f\n", RMSinduced);
                fprintf("NRMS value based on range: %f\n", NRMS_W_induced);
                fprintf("NRMS value based on absolute maximum: %f\n", NRMS_M_induced);
        
                % assigning for further analysis
                input_estimation = input_induced;
                output_estimation = output_induced;
                
                % Evaluation of validation data
                [estimation, deviation] = gp(hyp_opt, infv, meanfunc, covInduced, likfunc, input_estimation, output_estimation, input_validation); % extract the mean and covarance functions
                %% KPI-s
                [RMS, NRMS_W, NRMS_M, RMS_DEV] = KPIcalculation(estimation, deviation, output_validation);
                fprintf("EVALUTATION of snippet %d:\n", i);
                fprintf("RMS value is: %f\n", RMS);
                fprintf("NRMS value based on range: %f\n", NRMS_W);
                fprintf("NRMS value based on absolute maximum: %f\n", NRMS_M);
        
                % Evaluation of estimation data
                [estimationEstimation, deviationEstimation] = gp(hyp_opt, infv, meanfunc, covInduced, likfunc, input_estimation, output_estimation, input_estimation); % extract the mean and covarance functions
        
                %% KPI-s
                [RMSest, NRMS_W_est, NRMS_M_est, RMS_DEV_est] = KPIcalculation(estimationEstimation, deviationEstimation, output_estimation);
                fprintf("EVALUTATION of snippet %d:\n", i);
                fprintf("RMS_est value is: %f\n", RMSest);
                fprintf("NRMS_est value based on range: %f\n", NRMS_W_est);
                fprintf("NRMS_est value based on absolute maximum: %f\n", NRMS_M_est);
        
                % time series resimulation
                inputTimeSeries = (inputRaw-cin)./s_in;
                outputTimeSeries = (outputRaw-cout)./s_out(shiftID);
                [estimationTimeSeries, deviationTimeSeries] = gp(hyp_opt, infv, meanfunc, covInduced, likfunc, input_estimation, output_estimation, inputTimeSeries); % extract the mean and covarance functions
                offsets{1}.offset = estimationTimeSeries*s_out(1)+cout; offsets{1}.name = 'GPM'; offsets{1}.X = segment_m(:,indexes.X_abs); offsets{1}.marker = 'r';
                offsets{2}.offset = outputRaw; offsets{2}.name = 'REAL'; offsets{2}.X = segment_m(:,indexes.X_abs); offsets{2}.marker = 'k';
        
                f = plot_offsetPlots(segment_m,indexes,offsets, driverID);
                savefig(f, fullfile("plots",...
                    strcat('TimeSeriesPlot_shiftID_', num2str(shiftID), '_driver', num2str(driverID), '.fig')));
                saveas(f, fullfile("plots",...
                    strcat('TimeSeriesPlot_shiftID_', num2str(shiftID), '_driver', num2str(driverID), '.png')));
                close(f);

                plotDeviationPlots(estimation, deviation, output_validation, estimationEstimation, deviationEstimation, output_estimation, shiftID, driverID);
                
                trun = toc;
                KPI{shiftID}= [RMS NRMS_W NRMS_M RMS_DEV RMSest NRMS_W_est NRMS_M_est RMS_DEV_est hyp_opt.cov RMSinduced NRMS_W_induced NRMS_M_induced trun];
                ETA(shiftID).hyp_opt = hyp_opt;
                ETA(shiftID).input_estimation = input_estimation;
                ETA(shiftID).output_estimation = output_estimation;
                ETA(shiftID).normFactors.c_in = cin;
                ETA(shiftID).normFactors.s_in = s_in;
                ETA(shiftID).normFactors.c_out = cout;
                ETA(shiftID).normFactors.s_out = s_out;
                
                fprintf("time of shift id %d is %f\n", shiftID, trun);        
            end    
            save( fullfile("kpis",strcat('KPI_indSize_', num2str(indSize), '_driver_', num2str(driverID), '.mat')), 'KPI');
            save( fullfile("kpis",strcat('ETA_indSize_', num2str(indSize), '_driver_', num2str(driverID), '.mat')), 'ETA');
        end
    end
end

function plotDeviationPlots(estimation, deviation, output_validation, estimationEstimation, deviationEstimation, output_estimation, shiftID, driverID)
    f = figure('units','normalized','outerposition',[0 0 1 1]);
        
    % plot the estimation data
    subplot(2,1,1);
    confidenceBounds = [estimation+2*sqrt(deviation); flip(estimation-2*sqrt(deviation),1)];
    confidencePoints = (1:1:numel(estimation))';
    
    fill([confidencePoints; flip(confidencePoints)], confidenceBounds, 'y', 'DisplayName', '95% confidence');
    hold on;
    plot(estimation,'LineWidth', 2, 'color', 'k', 'LineStyle','--', 'DisplayName', 'estimated by GP');    
    ylabel('offset');
    grid on;
    plot(output_validation, 'LineWidth', 2, 'LineStyle', ':', 'color', 'b', 'DisplayName','validation data');
    legend;
    ylabel('offset(m)');
    
    % plot the estimation data
    subplot(2,1,2);
    confidenceBounds = [estimationEstimation+2*sqrt(deviationEstimation); flip(estimationEstimation-2*sqrt(deviationEstimation),1)];
    confidencePoints = (1:1:numel(estimationEstimation))';
    
    fill([confidencePoints; flip(confidencePoints)], confidenceBounds, 'y', 'DisplayName', '95% confidence');
    hold on;
    plot(estimationEstimation,'LineWidth', 2, 'color', 'k', 'LineStyle','--', 'DisplayName', 'estimated by GP');    
    ylabel('offset');
    grid on;
    plot(output_estimation, 'LineWidth', 2, 'LineStyle', ':', 'color', 'b', 'DisplayName','estimation data');
    legend;
    ylabel('offset(m)');
    
    savefig(f, fullfile("plots",...
                        strcat('ShiftID_', num2str(shiftID), '_driver', num2str(driverID), '.fig')));
    saveas(f, fullfile("plots",...
                    strcat('ShiftID_', num2str(shiftID), '_driver', num2str(driverID), '.png')));
    close(f);

end

function hyp = generateHyperParamStruct(KERNEL_TYPE, M)
    % ID problem
    if (KERNEL_TYPE == "@covLINiso")
        hyp = struct('mean', [], 'cov', 0, 'lik', -1); % initalize hyper-parameter structure associated with the mean, covariance an likelihood functions
    elseif (KERNEL_TYPE == "@covRQiso")
        hyp = struct('mean', [], 'cov', [0,0,0], 'lik', -1);
    elseif (KERNEL_TYPE=="{'covPERiso',{@covRQiso}}")
        hyp = struct('mean', [], 'cov', [0,0,0,0], 'lik', -1);
    elseif (KERNEL_TYPE=="{'covPERiso',{@covRQiso}}")
        hyp = struct('mean', [], 'cov', [0,0,0,0], 'lik', -1);
    elseif (KERNEL_TYPE == "{'covSum', {'covRQiso',{'covPERiso',{@covRQiso}}}}")
        hyp = struct('mean', [], 'cov', [0, 0, 0, 0, 0, 0, 0], 'lik', -1);
    elseif (KERNEL_TYPE == "{'covSum', {'covRQiso',{'covPERiso',{@covRQiso}}, 'covLINiso'}}")
        hyp = struct('mean', [], 'cov', [0, 0, 0, 0, 0, 0, 0, 0], 'lik', -1);
    elseif (KERNEL_TYPE == "{'covPERard', {@covSEard}}")
        hyp = struct('mean', [], 'cov', ones(1,3*M+1), 'lik', -1);
    elseif (KERNEL_TYPE == "{'covPPard',3}")
        hyp = struct('mean', [], 'cov', ones(1,M+1), 'lik', -1);
    elseif (KERNEL_TYPE == "{'covSum', {'covLINard', 'covSEard', {'covPPard',3}}}") 
        hyp = struct('mean', [], 'cov', ones(1,3*M+2), 'lik', -1);
    elseif (KERNEL_TYPE == "{'covSum', {'covSEiso', 'covLINiso'}}")
        hyp = struct('mean', [], 'cov', [0, 0, 0], 'lik', -1);
    elseif (KERNEL_TYPE == "{'covSum', {'covSEiso', 'covLINiso', {'covPERiso',{@covRQiso}}}}")
        hyp = struct('mean', [], 'cov', [0, 0, 0, 0, 0, 0, 0], 'lik', -1);
    elseif (KERNEL_TYPE == "{'covSum', {'covRQiso',{'covPERiso',{@covRQiso}}, {'covProd', {'covLINiso', 'covLINiso'}}}}")
        hyp = struct('mean', [], 'cov', [0, 0, 0, 0, 0, 0, 0, 0, 0], 'lik', -1);
    elseif (KERNEL_TYPE == "@covSEard" || KERNEL_TYPE == "{'covSEard'}")
        hyp = struct('mean', [], 'cov', ones(1,M+1), 'lik', -1);
    elseif (KERNEL_TYPE == "@covLINard")
        hyp = struct('mean', [], 'cov', ones(1,M), 'lik', -1);
    elseif (KERNEL_TYPE == "{'covSum', {'covSEard', 'covLINard'}}")
        hyp = struct('mean', [], 'cov', ones(1,2*M+1), 'lik', -1);
    elseif (KERNEL_TYPE == "{'covProd', {'covLINard', 'covLINard'}}")
        hyp = struct('mean', [], 'cov', ones(1,2*M), 'lik', -1);
    elseif (KERNEL_TYPE == "{'covSum',{{'covProd', {'covLINard', 'covLINard'}}, 'covSEard'}}")
        hyp = struct('mean', [], 'cov', ones(1,3*M+1), 'lik', -1);
    elseif (KERNEL_TYPE == "{'covSum', {'covLINard', {'covPPard',3}}}")
        hyp = struct('mean', [], 'cov', log(ones(1,2*M+1)), 'lik', -1);
    elseif (KERNEL_TYPE == "{'covSum', {'covSEard', {'covPPard',3}}}")
        hyp = struct('mean', [], 'cov', ones(1,2*M+2), 'lik', -1);
    else
        hyp = struct('mean', [], 'cov', [0, 0], 'lik', -1);
    end
end

function [RMS, NRMS_W, NRMS_M, RMS_DEV] = KPIcalculation(estimation, deviation, output)
    % RMS calculation
    RMS = sqrt(sum((estimation-output).^2)/size(estimation,1));
    % NRMS - normalization on scale
    W = max(output) - min(output);
    NRMS_W = RMS/W;
    % NRMS - normalization on absolute maximum
    M = max(abs(output));
    NRMS_M = RMS/M;
    % mean variance
    RMS_DEV = sqrt(sum(deviation.^2)/size(deviation,1));
end

function plotCorrelation (input, output, shiftID, driverID, variables)
    f = figure();
    set(f,'defaulttextInterpreter','latex') ;
    set(f, 'defaultAxesTickLabelInterpreter','latex');  
    set(f, 'defaultLegendInterpreter','latex');
    f.Position = [10 10 1000 1000];
    for i=1:size(input,2)
        subplot(4,2,i);
        plot(input(:,i), output, 'k.', 'MarkerSize', 4);
        grid on;
        xlim([-3,3]);
        ylim([-3,3]);
        title(strcat("$", variables(i),"\propto \delta$"));
        ylabel("$\delta$");
        xlabel(strcat("$", variables(i), "$"));
        set(gca,'FontSize', 14);
    end

    savefig(f, fullfile("plots",...
                        strcat('Correlation_plots_', num2str(shiftID), '_driver', num2str(driverID), '.fig')));
    saveas(f, fullfile("plots",...
                    strcat('Correlation_plots_', num2str(shiftID), '_driver', num2str(driverID), '.png')));
    close(f);
end
