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
DRIVER_ID = []; %[]; % [1, 5, 6, 10, 11, 13]; % possible values: 1,2,3,4,5...etc or []=all drivers simulated


%% DATA PIPELINE PARAMETERS
PARAMS.MAXIMUM_INPUT_SIZE = inf;
PARAMS.OUTPUT_SHIFT = linspace(15,150,10);
PARAMS.npDistances = [10, 39, 136];
p = 0.7;
dT = 0.1; % in seconds
K = floor(dT/0.05);
kappa_min = 2.5e-4;

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
    
    [~, ~, ~, ~, ~, s_out, ~, s_in] = prepareDataELDM(segmentMerged_m, indexesMerged, PARAMS);
    
    set(0,'DefaultFigureVisible','off');
        
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

        % removing some samples
        segment_m = segment_m(1:K:end,:);
        
        [~, ~, input, output, cout, ~, cin, ~] = prepareDataELDM(segment_m, indexes, PARAMS);
        
        inputRaw = input; % saved for further analysis
        outputRaw = output; % saved for further analysis

        % LIMIT DATA IF NEEDED
        % this is done before norm and central
        input = input(1:min(size(input,1),PARAMS.MAXIMUM_INPUT_SIZE),:);
        output = output(1:min(size(input,1),PARAMS.MAXIMUM_INPUT_SIZE), :);

        % NORM AND CENTRAL
        % According to the global cin and sin, cout and sout
        output = (output-cout)./s_out;
        input = (input-cin)./s_in;

        % SHUFFLE
        N = size(input,1);
        shuffledIndeces = randperm(N);
        input = input(shuffledIndeces,:);
        output = output(shuffledIndeces, :);
        
        % REMOVING SPARE POINTS I.E., OUTPUT HIGHER THAN 2, = 2 sigma
        % variance
        input = input(abs(output(:,1))<=3,:);
        output = output(abs(output(:,1))<=3,:);

        % DIVIDING BETWEEN TRAIN AND VALIDATION
        N = size(input,1);
        estimationData = 1:1:round(p*N);
        validationData = round(p*N)+1:1:N;
        input_train = input(estimationData,:,:);
        output_train = output(estimationData,:,:);
        input_validation = input(validationData,:,:);
        output_validation = output(validationData,:,:);

        tic;

        %% LDM model
        P = output_train'*input_train*inv(input_train'*input_train);
        trun_LDM = toc;
    
        y_train = input_train*P';
        y_val = input_validation*P';

        %% ELDM model
        % filter for left/right curves, no need to filter straight line, as
        % the output is centralized!
        tic;
        leftCurves = (mean((input_train.*s_in)') > kappa_min)';
        rightCurves = (mean((input_train.*s_in)') < -kappa_min)';

        inputTrainLeft = input_train(leftCurves,:);
        outputTrainLeft = output_train(leftCurves,:);
        inputTrainRight = input_train(rightCurves,:);
        outputTrainRight = output_train(rightCurves,:);

        P_left = outputTrainLeft'*inputTrainLeft*inv(inputTrainLeft'*inputTrainLeft);
        P_right = outputTrainRight'*inputTrainRight*inv(inputTrainRight'*inputTrainRight);

        y_ELDM_train = zeros(size(input_train,1),3);
        y_ELDM_train(leftCurves,:) = inputTrainLeft*P_left';
        y_ELDM_train(rightCurves,:) = inputTrainRight*P_right';

        leftCurves = mean(input_validation) > kappa_min;
        rightCurves = mean(input_validation) < -kappa_min;

        inputValidationLeft = input_validation(leftCurves,:);
        outputValidationLeft = input_validation(leftCurves,:);
        inputValidationRight = input_validation(rightCurves,:);
        outputValidationRight = input_validation(rightCurves,:);

        y_ELDM_val = zeros(size(input_validation,1),3);
        y_ELDM_val(leftCurves,:) = inputValidationLeft*P_left';
        y_ELDM_val(rightCurves,:) = inputValidationRight*P_left';

        trun_ELDM = toc;

        for shiftID=1:size(output,2)
            % LDM KPIs
            [RMS_LDM, NRMS_W_LDM, NRMS_M_LDM] = KPIcalculation(y_train(:,shiftID), output_train(:,shiftID));
            [RMS_LDM_val, NRMS_LDM_W_val, NRMS_LDM_M_val] = KPIcalculation(y_val(:,shiftID), output_validation(:,shiftID));    

            % ELDM KPIs
            [RMS_ELDM, NRMS_W_ELDM, NRMS_M_ELDM] = KPIcalculation(y_ELDM_train(:,shiftID), output_train(:,shiftID));
            [RMS_ELDM_val, NRMS_W_ELDM_val, NRMS_M_ELDM_val] = KPIcalculation(y_ELDM_val(:,shiftID), output_validation(:,shiftID));
    
            KPI{shiftID} = [RMS_LDM, NRMS_W_LDM, NRMS_M_LDM, 0, RMS_LDM_val, NRMS_LDM_W_val, NRMS_LDM_M_val, 0, trun_LDM, ...
                            RMS_ELDM, NRMS_W_ELDM, NRMS_M_ELDM, 0, RMS_ELDM_val, NRMS_W_ELDM_val, NRMS_M_ELDM_val, 0, trun_ELDM];
        end

        ETA.P = P;
        ETA.P_left = P_left;
        ETA.P_right = P_right;
        ETA.input_validation = input_validation;
        ETA.output_validation = output_validation;
        ETA.input_train = input_train;
        ETA.output_train = output_train;
        ETA.input_trainLeft = inputTrainLeft;
        ETA.input_trainRight = inputTrainRight;
        ETA.input_validationLeft = inputValidationLeft;
        ETA.input_validationRight = inputValidationRight;
        ETA.normFactors.c_in = cin;
        ETA.normFactors.s_in = s_in;
        ETA.normFactors.c_out = cout;
        ETA.normFactors.s_out = s_out;
 
        save( fullfile("kpis",strcat('KPI_driver_', num2str(driverID), '.mat')), 'KPI');
        save( fullfile("kpis",strcat('ETA_driver_', num2str(driverID), '.mat')), 'ETA');

        clear KPI
    end
end

function [RMS, NRMS_W, NRMS_M] = KPIcalculation(estimation, output)
    % RMS calculation
    RMS = sqrt(sum(sum((estimation-output).^2))/(size(estimation,1)*size(estimation,2)));
    % NRMS - normalization on scale
    W = max(max(output)) - min(min(output));
    NRMS_W = RMS/W;
    % NRMS - normalization on absolute maximum
    M = max(max(abs(output)));
    NRMS_M = RMS/M;
end