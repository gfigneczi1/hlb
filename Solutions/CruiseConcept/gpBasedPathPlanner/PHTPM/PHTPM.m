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
DRIVER_ID = 1; % [1, 5, 6, 10, 11, 13]; % possible values: 1,2,3,4,5...etc or []=all drivers simulated
ITERATE_INDUCTION = true;

%% DATA PIPELINE PARAMETERS
MAXIMUM_INPUT_SIZE = inf; % before snipetting and normalization
KERNEL_TYPE = "{'covPPard',3}";
RATIO_OF_TRAIN_VS_TOTAL = 0.7;
MAXIMUM_PREVIEW_DISTANCE = 150; % from within preview information is extracted, unit: m
OUTPUT_SHIFT = 0; %possible values: scalar number or array e.g., linspace(15,MAXIMUM_PREVIEW_DISTANCE,10);
INDUCED_SIZE = 100; % only applicable if ITERATE_INDUCTION is FALSE
Lp = 1.5; % in seconds
HISTORICAL_WINDOW_LENGTH = 15; % for 0.1 s sampling time based on paper 
PREDICTION_WINDOW_LENGTH = 1.5; % preview horizon in seconds

variables = ["t_{ot}", "f_{ot}", "v_x", "a_x", "\omega", "\kappa", "d\kappa"];

if (isempty(DRIVER_ID))
    i0 = 1; i1 = length(segments.segments);
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
    PARAMS.Lp = Lp;
    [~, ~, ~, ~, ~, s_out, ~, s_in] = prepareDataPHTPM(segmentMerged_m, indexesMerged, PARAMS);
    
    set(0,'DefaultFigureVisible','off');
    
    % PARAMETERS REMAPPING
    shiftOnOutputSelection = OUTPUT_SHIFT;  %shift the offset in time (positive means shift forward)
    p = RATIO_OF_TRAIN_VS_TOTAL; %percentage of evaluation data from entire dataset
    
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

        % resampling to have dT = 0.1 s
        if (dT < 0.1)
            k = round(0.1/dT);
        end
        segment_m = segment_m (1:k:end,:);

        fp = PREDICTION_WINDOW_LENGTH/0.1;

        tic;
        % calculate thetaTP and thetaFP
        ladVector = Lp*segment_m(:,indexes.VelocityX);
        ye = segment_m(:,indexes.c0)+ ...
            tan(segment_m(:,indexes.LaneOrientation)).*ladVector + ...
            segment_m(:,indexes.LaneCurvature)/2.*ladVector.^2 + ...
            segment_m(:,indexes.c3)/6.*ladVector.^3;
        thetaFP = atan(ye/Lp);

        determinant = 4*(segment_m(:,indexes.LaneCurvature)/2).^2 - 12*tan(segment_m(:,indexes.LaneOrientation)).*segment_m(:,indexes.c3)/6;
        determinant(determinant<0) = nan;

        x1 = (-2*segment_m(:,indexes.LaneCurvature)/2 + determinant)./(segment_m(:,indexes.c3));
        x2 = (-2*segment_m(:,indexes.LaneCurvature)/2 - determinant)./(segment_m(:,indexes.c3));
        x = min(max(0,x1),max(0,x2));
        x(isnan(determinant)) = Lp*segment_m(isnan(determinant), indexes.VelocityX);
        x(x==0) = Lp*segment_m(x==0, indexes.VelocityX);

        kappa_r = segment_m(:,indexes.LaneCurvature)+segment_m(:,indexes.c3).*x;
        thetaTP = atan(x.*kappa_r);

        input = [movmean(thetaTP,100), ...
                movmean(thetaFP,100), ...
                segment_m(:, indexes.VelocityX), ...
                movmean(segment_m(:, indexes.Acceleration_Y),100), ...
                movmean(segment_m(:, indexes.YawRate),100), ...
                movmean(-segment_m(:, indexes.LaneOrientation), 20)];
        % note: lane orientation is multiplied by -1 to yield the
        % right meaning of the phi value given in the paper
        
        inputRaw = input; % saved for further analysis

        % Laterial deviation index
        clear output
        output = -0.5*(segment_m(:,indexes.LaneEdgePositionLeft)+segment_m(:,indexes.LaneEdgePositionRight));

        outputRaw = output; % saved for further analysis

        % LIMIT DATA IF NEEDED
        % this is done before norm and central
        input = input(1:min(size(input,1),MAXIMUM_INPUT_SIZE),:);
        output = output(1:min(size(input,1),MAXIMUM_INPUT_SIZE), :);

        % NORM AND CENTRAL
        % According to the global cin and sin, cout and sout
        [~, cout, ~] = normalize(output(:,1));
        [~, cin, ~] = normalize(input);
        output(:,1) = (output(:,1)-cout)./s_out;
        input = (input-cin)./s_in;
        
        % REMOVING SPARE POINTS I.E., OUTPUT HIGHER THAN 2, = 2 sigma
        % variance
        input = input(abs(output(:,1))<=3,:);
        output = output(abs(output(:,1))<=3,:);

        %% Define PHTPM
        % producing the h-window of the input data
        for n = HISTORICAL_WINDOW_LENGTH+1:size(input,1)-fp+1
            if (n==HISTORICAL_WINDOW_LENGTH+1)
                X(n-HISTORICAL_WINDOW_LENGTH,:,:) = [input(1:n-1,:), output(1:n-1,1)]';
                T(n-HISTORICAL_WINDOW_LENGTH,:,:) = output(n:n+fp-1,1)';
            else
                X(n-HISTORICAL_WINDOW_LENGTH,:,:) = [input(n-HISTORICAL_WINDOW_LENGTH:n-1,:), output(n-HISTORICAL_WINDOW_LENGTH:n-1,1)]';
                T(n-HISTORICAL_WINDOW_LENGTH,:,:) = output(n:n+fp-1,1)';
            end
        end
        
        % SHUFFLE
        N = size(X,1);
        shuffledIndeces = randperm(N);
        X = X(shuffledIndeces,:,:);
        T = T(shuffledIndeces,:,:);

        % DIVIDING BETWEEN TRAIN AND VALIDATION
        M = size(X,1);
        % EVALUATION / VALIDATION DATA SELECTION
        estimationData = 1:1:round(p*M);
        validationData = round(p*M)+1:1:M;
        XTrain_ = X(estimationData,:,:);
        TTrain_ = T(estimationData,:,:);
        XVal_ = X(validationData,:,:);
        TVal_ = T(validationData,:,:);

        for n=1:size(XTrain_,1)
            XTrain{n,1}(:,:) = XTrain_(n,:,:);
            TTrain{n,1}(1,:) = TTrain_(n,:,:);
        end
        for n=1:size(XVal_,1)
            XVal{n,1}(:,:) = XVal_(n,:,:);
            TVal{n,1}(1,:) = TVal_(n,:,:);
        end

        numResponses = size(TTrain{1,1},1);
        numFeatures = size(XTrain{1,1},1);

        layers = [
            sequenceInputLayer(numFeatures);
            lstmLayer(100, OutputMode="last");
            fullyConnectedLayer(numFeatures);
            lstmLayer(50, OutputMode="sequence");
            fullyConnectedLayer(numResponses);
            regressionLayer
            ];

        
        TVal__(:,:) = TVal_(:,1,:);

        options = trainingOptions("adam", ...
        MaxEpochs=50, ...
        InitialLearnRate=0.05,...
        Shuffle="never", ...
        ValidationData = {XVal, TVal}, ...
        GradientThreshold=1, ...
        Verbose=false, ...
        ValidationFrequency = 1, ...
        L2Regularization = 0.0001, ...
        LearnRateDropPeriod=125, ...
        LearnRateDropFactor=0.1, ...
        MiniBatchSize = 620, ...
        Plots="training-progress");

        TTrain__(:,:) = TTrain_(:,1,:);
        net = trainNetwork(XTrain,TTrain,layers,options);
                
        % Evaluation of validation data
        YVal = predict(resetState(net), XVal);
        for n=1:size(YVal,1)
            Yval_(n,:) = YVal{n,1}(1,:);
            Yvalref_(n,:) = TVal{n,1}(1,:);
        end

        % Evaluation of estimation data
        YEst = predict(resetState(net), XTrain);
        for n=1:size(YEst,1)
            YEst_(n,:) = YEst{n,1}(1,:);
            YEstref_(n,:) = TTrain{n,1}(1,:);
        end

        trun = toc;
        %% KPI-s
        for shiftID = 1:size(Yvalref_,2)
            [RMS, NRMS_W, NRMS_M] = KPIcalculation(Yval_(:,shiftID), Yvalref_(:,shiftID));
            [RMSest, NRMS_W_est, NRMS_M_est] = KPIcalculation(YEst_(:,shiftID), YEstref_(:,shiftID));

            KPI{shiftID} = [RMS NRMS_W NRMS_M 0 RMSest NRMS_W_est NRMS_M_est 0 trun/size(Yvalref_,2)];
        end
        
        fprintf("EVALUTATION of driver %d:\n", driverID);
        fprintf("RMS value is: %f\n", RMS);
        fprintf("NRMS value based on range: %f\n", NRMS_W);
        fprintf("NRMS value based on absolute maximum: %f\n", NRMS_M);        

        %% KPI-s
        
        fprintf("EVALUTATION of driver %d:\n", driverID);
        fprintf("RMS_est value is: %f\n", RMSest);
        fprintf("NRMS_est value based on range: %f\n", NRMS_W_est);
        fprintf("NRMS_est value based on absolute maximum: %f\n", NRMS_M_est);
        
        ETA.net = net;
        ETA.XTrain = XTrain;
        ETA.TTrain = TTrain;
        ETA.XVal = XVal;
        ETA.TVal = TVal;
        ETA.normFactors.c_in = cin;
        ETA.normFactors.s_in = s_in;
        ETA.normFactors.c_out = cout;
        ETA.normFactors.s_out = s_out;
        
        fprintf("time of driver %d is %f\n", driverID, trun);        
        save( fullfile("kpis",strcat('KPI_driver_', num2str(driverID), '.mat')), 'KPI');
        save( fullfile("kpis",strcat('ETA_driver_', num2str(driverID), '.mat')), 'ETA');
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
                        strcat('ShiftID_', num2str, '_driver', num2str(driverID), '.fig')));
    saveas(f, fullfile("plots",...
                    strcat('ShiftID_', num2str, '_driver', num2str(driverID), '.png')));
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

function [RMS, NRMS_W, NRMS_M] = KPIcalculation(estimation, output)
    % RMS calculation
    RMS = sqrt((sum((estimation-output).^2))/(size(estimation,1)));
    % NRMS - normalization on scale
    W = (max(output)) - (min(output));
    NRMS_W = RMS/W;
    % NRMS - normalization on absolute maximum
    M = (max(abs(output)));
    NRMS_M = RMS/M;
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
                        strcat('Correlation_plots_', num2str, '_driver', num2str(driverID), '.fig')));
    saveas(f, fullfile("plots",...
                    strcat('Correlation_plots_', num2str, '_driver', num2str(driverID), '.png')));
    close(f);
end
