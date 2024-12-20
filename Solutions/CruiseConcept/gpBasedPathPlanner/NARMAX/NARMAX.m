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


%% DATA PIPELINE PARAMETERS
PARAMS.MAXIMUM_INPUT_SIZE = inf;
PARAMS.OUTPUT_SHIFT = linspace(15,150,10);
p = 0.7;
dT = 0.1; % in seconds
K = floor(dT/0.05);

variables = ["t_{ot}", "f_{ot}", "v_x", "a_x", "\omega", "\kappa", "d\kappa"];

if (isempty(DRIVER_ID))
    i0 = 12; i1 = length(segments.segments);
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
    
    [~, ~, ~, ~, ~, s_out, ~, s_in] = prepareDataNARMAX(segmentMerged_m, indexesMerged, PARAMS);
    
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
        
        [~, ~, input, output, cout, ~, cin, ~] = prepareDataNARMAX(segment_m, indexes, PARAMS);
        
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

        %% ARX model
        na = 2; nb = ones(1,7)*2; nk = ones(1,7);
        sysARX = arx(input_train, output_train(:,1),[na nb nk],InputName={'u1', 'u2', 'u3', 'u4', 'u5', 'u6', 'u7'}, OutputName={'y'}, Ts=dT);

        y = predict(sysARX, [output_train(:,1) input_train], inf);
        [RMS_ARX, NRMS_W_ARX, NRMS_M_ARX] = KPIcalculation(y, output_train(:,1));

        y = predict(sysARX, [output_validation(:,1) input_validation], inf);
        [RMS_ARX_val, NRMS_W_ARX_val, NRMS_M_ARX_val] = KPIcalculation(y, output_validation(:,1));

        trun_ARX = toc;
        tic;

        %% NARX model
        G = idGaussianProcess('SquaredExponential');
        S = idSigmoidNetwork;
        opt = nlarxOptions(Focus="simulation", Display="on",Normalize=false);
        opt.SearchOptions.MaxIterations = 50;

        output_lag = [1,2];
        input_lag = [1,2,3];

        names = ["y", "u1", "u2", "u3", "u4", "u5", "u6", "u7"];
        lags = {output_lag, input_lag, input_lag, input_lag, input_lag, input_lag, input_lag, input_lag};
        lreg = linearRegressor(names,lags);
        %sysNARX_sigmoid = nlarx([output_train(:,1) input_train], lreg, S, opt)
        sysNARX_gp = nlarx([output_train(:,1) input_train], lreg, G, opt);

        y = predict(sysNARX_gp, [output_train(:,1) input_train], inf);
        [RMS_NARX, NRMS_W_NARX, NRMS_M_NARX] = KPIcalculation(y, output_train(:,1));

        y = predict(sysNARX_gp, [output_validation(:,1) input_validation], inf);
        [RMS_NARX_val, NRMS_W_NARX_val, NRMS_M_NARX_val] = KPIcalculation(y, output_validation(:,1));

        trun_NARX = toc;
        trun = toc;
        KPI = [RMS_ARX, NRMS_W_ARX, NRMS_M_ARX, 0, RMS_ARX_val, NRMS_W_ARX_val, NRMS_M_ARX_val, trun_ARX; ...
                RMS_NARX, NRMS_W_NARX, NRMS_M_NARX, 0, RMS_NARX_val, NRMS_W_NARX_val, NRMS_M_NARX_val, trun_NARX];
        ETA.sysARX = sysARX;
        ETA.sysNARX = sysNARX_gp;
        ETA.input_validation = input_validation;
        ETA.output_validation = output_validation;
        ETA.input_train = input_train;
        ETA.output_train = output_train;
        ETA.normFactors.c_in = cin;
        ETA.normFactors.s_in = s_in;
        ETA.normFactors.c_out = cout;
        ETA.normFactors.s_out = s_out;

        KPI_sum{driverID} = KPI;
        ETA_sum{driverID} = ETA;

        save( fullfile("kpis",strcat('KPI_driver_', num2str(driverID), '.mat')), 'KPI');
        save( fullfile("kpis",strcat('ETA_driver_', num2str(driverID), '.mat')), 'ETA');
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

function y = ARX(yk, uk, ek, A, B, C)
    % A(z)*y(k) = B(z)*u(k)+C(z)*e(k)
    % rearrange: y(k) = B(z)*u(k)+C(z)*e(k)-A(z-1)y(k-1)
    % numel(A) == na, ...etc
    y = B*uk + C*ek - A*yk;
end
