close all; clc; clear;
load('C:\database\KDP_HLB_GP\Dr008_Dr023_input_GP_fixedAcceleration.mat');

addpath(fullfile("..", "library"));
addpath(fullfile("..", "library", "gpml-matlab-master"));
addpath(fullfile("..", "library", "gpml-matlab-master", "cov"));
addpath(fullfile("..", "library", "gpml-matlab-master", "inf"));
addpath(fullfile("..", "library", "gpml-matlab-master", "lik"));
addpath(fullfile("..", "library", "gpml-matlab-master", "mean"));
addpath(fullfile("..", "library", "gpml-matlab-master", "prior"));
addpath(fullfile("..", "library", "gpml-matlab-master", "util"));

PARAMS.MAXIMUM_INPUT_SIZE = 15000; % before snipetting and normalization
PARAMS.KERNEL_TYPE ="{'covPPard',3}"; % "{'covSum', {'covLINard', {'covPPard',3}}}"; %"{'covSum', {'covSEard', {'covPPard',3}}}";
PARAMS.RATIO_OF_TRAIN_VS_TOTAL = 0.7;
PARAMS.MAXIMUM_PREVIEW_DISTANCE = 150; % from within preview information is extracted, unit: m
PARAMS.OUTPUT_STEP_DISTANCE = 10; % in meters
PARAMS.NUMBER_OF_PREVIEW_INFORMATION = 2; % maximum number
PARAMS.OUTPUT_SHIFT = linspace(15,PARAMS.MAXIMUM_PREVIEW_DISTANCE,10); %10:OUTPUT_STEP_DISTANCE:MAXIMUM_PREVIEW_DISTANCE; % preview distance is divided into sections where the output will be predicted
PARAMS.EPOCH_CALCULATION = true;
PARAMS.EPOCH_NUMBER = 10;
PARAMS.GREEDY_REDUCTION = true;
PARAMS.LDM_NP = [10, 39, 136];
PARAMS.DriverID = 4;
PARAMS.InductionID = 4;
PARAMS.FILTER_OUTPUT = false;

for driverID=1:8
    PARAMS.DriverID = driverID;
    segment = segments.segments(driverID).segment;
    name = segments.segments(driverID).name;
    [~, segment_m, indexes] = prepareInputForPlanner(segment);
    [estimation, deviation, inputRaw, outputRaw, input, output] = gpGenerateEstimate(segment_m, indexes, PARAMS);

    for j=1:size(estimation,2)
        plotDeviationPlots(estimation(:,j), deviation(:,j), outputRaw(:,j), j, driverID, segment_m(:,indexes.X_abs));
    end
    
end

function [estimation, deviation, inputRaw, outputRaw, input, output] = gpGenerateEstimate(segment_m, indexes, PARAMS)
    %% step 0: generate input data
    [~, ~, inputRaw, outputRaw] = prepareData(segment_m, indexes, PARAMS);
    
    %% step 1: parameter loading
    % read the learnt data from previous runs
    pathToParams = "C:\git\KDP\publications\GP\results\sparseGP_fixedAcceleration\kpis";
    paramGP = dir(fullfile(pathToParams,"ETA_*"));

    for fileID = 1:length(paramGP)
        paramDriverID = str2num(paramGP(fileID).name(strfind(paramGP(fileID).name, 'driver_')+7:strfind(paramGP(fileID).name, '.')-1));
        inductionSizeID = str2num(paramGP(fileID).name(strfind(paramGP(fileID).name, 'indSize_')+8:strfind(paramGP(fileID).name, 'indSize_')+8));
        if (paramDriverID == PARAMS.DriverID && inductionSizeID == PARAMS.InductionID)
            paramData = load(fullfile(paramGP(fileID).folder, paramGP(fileID).name));
            for npID=1:10
                GP_params.hyp_opt_array{npID} = paramData.ETA(npID).hyp_opt;
                GP_params.output_estimation{npID} = paramData.ETA(npID).output_estimation;
                GP_params.input_estimation{npID} = paramData.ETA(npID).input_estimation;
            end
            if (isfield(paramData.ETA(npID).normFactors, "c_in"))
                c_in = paramData.ETA(npID).normFactors.c_in;
                s_in = paramData.ETA(npID).normFactors.s_in;
                c_out(1:10) = paramData.ETA(npID).normFactors.c_out;
                s_out = paramData.ETA(npID).normFactors.s_out;
            else
                c_in = paramData.ETA(npID).normFactors(1:7); % common for all GPs
                s_in = paramData.ETA(npID).normFactors(8:14); % common for all GPs
                c_out(1:10) = paramData.ETA(npID).normFactors(15); % one at each node point
                s_out = paramData.ETA(npID).normFactors(16:25); % one at each node point
            end
            break;
        else
            paramData = [];
        end
    end

    %% step 3: norm and central
    for i=1:size(inputRaw,2)
        input(:,i) = (inputRaw(:,i)-c_in(i))/s_in(i);
    end
    for i=1:size(outputRaw,2)
        output(:,i) = (outputRaw(:,i)-c_out(i))/s_out(i);
    end

    %% step 4: generate output
    meanfunc = [];       % Start with a zero mean prior
    eval(strcat('covfunc = ',PARAMS.KERNEL_TYPE));    % Squared Exponental covariance function
    % ID problem            
    likfunc = @likGauss;    % Gaussian likelihood
    for i = 1:length(GP_params.hyp_opt_array)
        hyp = struct('mean', [], 'cov', 0, 'lik', -1);
        hyp.cov = GP_params.hyp_opt_array{i}.cov;
        hyp.lik = GP_params.hyp_opt_array{i}.lik;
        [estimationGP(:,i), deviationGP(:,i)] = gp(hyp, @infGaussLik, meanfunc, covfunc, likfunc, GP_params.input_estimation{i}, GP_params.output_estimation{i}, input); % extract the mean and covarance functions
    end

    %% step 5: node point calculation
    for i=1:size(estimationGP,2)
        estimation(:,i) = estimationGP(:,i)*s_out(i)+c_out(i);
        deviation(:,i) = 2*sqrt(deviationGP(:,i))*s_out(i);
    end
end

function plotDeviationPlots(estimation, deviation, output_validation, shiftID, driverID, X)
    f = figure('units','normalized','outerposition',[0 0 1 1]);
        
    % plot the estimation data
    confidenceBounds = [estimation+deviation; flip(estimation-deviation,1)];
    confidencePoints = X; %(1:1:numel(estimation))';
    
    fill([confidencePoints; flip(confidencePoints)], confidenceBounds, 'y', 'DisplayName', '95% confidence');
    hold on;
    plot(X, estimation,'LineWidth', 2, 'color', 'k', 'LineStyle','--', 'DisplayName', 'estimated by GP');    
    ylabel('offset');
    grid on;
    plot(X, output_validation, 'LineWidth', 2, 'LineStyle', ':', 'color', 'b', 'DisplayName','validation data');
    legend;
    ylabel('offset(m)');  

    savefig(f, fullfile("plots",...
                        strcat('ShiftID_', num2str(shiftID), '_driver', num2str(driverID), '.fig')));
    saveas(f, fullfile("plots",...
                    strcat('ShiftID_', num2str(shiftID), '_driver', num2str(driverID), '.png')));
    close(f);
end