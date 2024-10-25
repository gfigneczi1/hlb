% This function is the cover function to call the GP sub-function for a
% given segment_m. segment_ms is an array of struct containing the driving data
% of multiple drivers. Config is a metadata struct containing information
% about e.g., the plot folders
close all; clear;
load('C:\database\KDP_HLB_GP\Dr008_Dr023_input_GP.mat');
config.root = "./";

PARAMS.MAXIMUM_INPUT_SIZE = 15000; % before snipetting and normalization
PARAMS.KERNEL_TYPE = "{'covPPard',3}";
PARAMS.RATIO_OF_TRAIN_VS_TOTAL = 0.5;
PARAMS.MAXIMUM_PREVIEW_DISTANCE = 150; % from within preview information is extracted, unit: m
PARAMS.OUTPUT_STEP_DISTANCE = 10; % in meters
PARAMS.NUMBER_OF_PREVIEW_INFORMATION = 2; % maximum number
PARAMS.OUTPUT_SHIFT = linspace(15,PARAMS.MAXIMUM_PREVIEW_DISTANCE,10); %10:OUTPUT_STEP_DISTANCE:MAXIMUM_PREVIEW_DISTANCE; % preview distance is divided into sections where the output will be predicted
PARAMS.FILTER_OUTPUT = false;

temp_folder_path = config.root;
plots_folder_name = 'plots';
set(0,'DefaultFigureVisible','off');

% PARAMETERS
shiftOnOutputSelection = PARAMS.OUTPUT_SHIFT;  %shift the offset in time (positive means shift forward)
p = PARAMS.RATIO_OF_TRAIN_VS_TOTAL; %percentage of evaluation data from entire dataset
PARAMS.DRIVERID_GROUP = [1, 6];

stochasticClusters = clusterDrivers();

for refDriverID = 1:2
    pathToParams = "C:\git\KDP\publications\GP\results\sparseGP_chosenParams\kpis";
    paramGP = dir(fullfile(pathToParams,"ETA_*"));
    driverID_ = PARAMS.DRIVERID_GROUP(refDriverID);
    for fileID = 1:length(paramGP)
        paramDriverID = str2num(paramGP(fileID).name(strfind(paramGP(fileID).name, 'driver_')+7:strfind(paramGP(fileID).name, '.')-1));
        if (paramDriverID == driverID_)
            paramData = load(fullfile(paramGP(fileID).folder, paramGP(fileID).name));
            for npID=1:10
                GP_params.hyp_opt_array{npID} = paramData.ETA(npID).hyp_opt;
                GP_params.output_estimation{npID} = paramData.ETA(npID).output_estimation;
                GP_params.input_estimation{npID} = paramData.ETA(npID).input_estimation;
            end
            c_in = paramData.ETA(npID).normFactors(1:7); % common for all GPs
            c_out(1:10) = paramData.ETA(npID).normFactors(15); % one at each node point
            s_in = paramData.ETA(npID).normFactors(8:14); % common for all GPs
            s_out = paramData.ETA(npID).normFactors(16:25); % one at each node point
            break;
        end
    end
                        
    for driverID = 1:size(segments.segments,2)
        DRIVER_ID_IF_NOT_MERGED = driverID;
        segment = segments.segments(DRIVER_ID_IF_NOT_MERGED).segment;
        name = segments.segments(DRIVER_ID_IF_NOT_MERGED).name;
        % transforming the struct to array for lower calculation time
        [~, segment_m, indexes] = prepareInputForPlanner(segment);    
            
        [~, ~, inputRaw, outputRaw] = prepareData(segment_m, indexes, PARAMS);
        clear input output 
        for i=1:size(inputRaw,2)
            input(:,i) = (inputRaw(:,i)-c_in(i))/s_in(i);
        end
        for i=1:size(outputRaw,2)
            output(:,i) = (outputRaw(:,i)-c_out(i))/s_out(i);
        end
        
        [estimationGP, deviationGP] = resimulateTimeSequence (input, GP_params, PARAMS);

        % error calculation
        e = [];
        for i=1:length(estimationGP)
            e = [e; estimationGP{i}-output(:,i)];
        end
        error(refDriverID,driverID) = norm(e,2);
    
        %f = plotOffsets(estimationGP, deviationGP, segment_m, indexes, 1, outputRaw, s_in, c_in, s_out, c_out, driverID);
    end
end

%% SUPPORT FUNCTIONS
function [estimationGP, deviationGP, estimationLRM, estimationLDM] = resimulateTimeSequence (input_validation, GP_params, PARAMS)
    % GP resimulate
    meanfunc = [];       % Start with a zero mean prior
    eval(strcat('covfunc = ',PARAMS.KERNEL_TYPE));    % Squared Exponental covariance function
    % ID problem            
    likfunc = @likGauss;    % Gaussian likelihood
    for i = 1:length(GP_params.hyp_opt_array)
        hyp = struct('mean', [], 'cov', 0, 'lik', -1);
        hyp.cov = GP_params.hyp_opt_array{i}.cov;
        hyp.lik = GP_params.hyp_opt_array{i}.lik;
        [estimationGP{i}, deviationGP{i}] = gp(hyp, @infGaussLik, meanfunc, covfunc, likfunc, GP_params.input_estimation{i}, GP_params.output_estimation{i}, input_validation); % extract the mean and covarance functions
    end
end

function c = clusterDrivers()
    pathToParams = "C:\git\KDP\publications\GP\results\sparseGP_chosenParams\kpis";
    paramGP = dir(fullfile(pathToParams,"ETA_*"));
    for fileID = 1:length(paramGP)
        paramDriverID = str2num(paramGP(fileID).name(strfind(paramGP(fileID).name, 'driver_')+7:strfind(paramGP(fileID).name, '.')-1));
        paramData = load(fullfile(paramGP(fileID).folder, paramGP(fileID).name));
        N = length(paramData.ETA(1).hyp_opt.cov)-1; % -1: sigma_f is excluded
        for npID=1:10
            GP_params((npID-1)*N+1:npID*N, paramDriverID) = exp(paramData.ETA(npID).hyp_opt.cov(1:end-1))'; % columns are the parameters of one driver
        end
    end
    % cluster
    rng(1);
    c = kmeans(GP_params',3, 'Distance','sqeuclidean', 'MaxIter', 1000, 'Replicates',10);
end


function f = plotOffsets(estimationGP, deviationGP, segment_m, indexes, i, outputRaw, s_in, c_in, s_out, c_out, driverID)

    relSegment = [1.4185e6, 1.4205e6];
    dx = mean(diff(segment_m(:,indexes.X_abs)));
    X_abs = segment_m(:,indexes.X_abs);

    f = figure(i);
    f.Position = [10 10 1000 1500];
    set(f,'defaulttextInterpreter','latex') ;
    set(f, 'defaultAxesTickLabelInterpreter','latex');  
    set(f, 'defaultLegendInterpreter','latex');

    subplot(3,1,1);
    straightSections = abs(segment_m(:,indexes.LaneCurvature))<2.5e-4;
    straightBoundingBoxes = zeros(numel(find(diff(straightSections) < 0))+numel(find(diff(straightSections) > 0))+1,1);
    if (straightSections(1))
        % start with straight
        straightBoundingBoxes(2:2:end) = find(diff(straightSections) < 0);
        straightBoundingBoxes(1:2:end) = [1; find(diff(straightSections) > 0)];
    elseif (straightSections(end))
        % end with straight
        straightBoundingBoxes(2:2:end) = [find(diff(straightSections) < 0); size(segment_m,1)];
        straightBoundingBoxes(1:2:end) = [find(diff(straightSections) > 0)];
    else
        straightBoundingBoxes(1:2:end-1) = find(diff(straightSections) > 0);
        straightBoundingBoxes(2:2:end-1) = find(diff(straightSections) < 0);
    end

    confidenceBounds = [estimationGP{i}+2*sqrt(deviationGP{i}); flip(estimationGP{i}-2*sqrt(deviationGP{i}),1)];
    fill([X_abs; flip(X_abs)], confidenceBounds*s_out(i)+c_out(i), 'y', 'DisplayName', '95\% confidence - GP');
    grid on; hold on;
    for j=1:length(straightBoundingBoxes)/2
        if (j==1)
            fill([X_abs(straightBoundingBoxes(2*j-1)); X_abs(straightBoundingBoxes(2*j)); X_abs(straightBoundingBoxes(2*j)); X_abs(straightBoundingBoxes(2*j-1))], ...
                [-1;-1;1;1], 'g', 'FaceAlpha', 0.2, 'DisplayName', 'Straight sections', 'LineStyle','none');
        else
            fill([X_abs(straightBoundingBoxes(2*j-1)); X_abs(straightBoundingBoxes(2*j)); X_abs(straightBoundingBoxes(2*j)); X_abs(straightBoundingBoxes(2*j-1))], ...
                [-1;-1;1;1], 'g', 'FaceAlpha', 0.2, 'LineStyle','none', "HandleVisibility","off");
        end
    end
    plot(X_abs, estimationGP{i}*s_out(i)+c_out(i), 'color', 'r',  'LineWidth', 2, 'DisplayName', 'GP');
    plot(X_abs, outputRaw(:,i), 'color', 'k', 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'Reference');
    xlabel('X-UTM(m)'); ylabel("$\delta(m)$");
    yline(0, "HandleVisibility","off", 'Alpha',0.3, 'Color','k', 'LineWidth',3);

    legend("Location", "best", "Orientation","horizontal");
    set(gca,'FontSize', 14);
    xlim(relSegment);
    ylim([-1,1]);
    title(strcat("Estimation of $\delta_1$ - Driver", {' '}, num2str(driverID)));

    subplot(3,1,2);
    plot(X_abs, segment_m(:,indexes.LaneCurvature), 'LineWidth', 2, 'color', 'k');
    grid on; hold on;
    for j=1:length(straightBoundingBoxes)/2
        if (j==1)
            fill([X_abs(straightBoundingBoxes(2*j-1)); X_abs(straightBoundingBoxes(2*j)); X_abs(straightBoundingBoxes(2*j)); X_abs(straightBoundingBoxes(2*j-1))], ...
                [-1;-1;1;1], 'g', 'FaceAlpha', 0.2, 'DisplayName', 'Straight sections', 'LineStyle','none');
        else
            fill([X_abs(straightBoundingBoxes(2*j-1)); X_abs(straightBoundingBoxes(2*j)); X_abs(straightBoundingBoxes(2*j)); X_abs(straightBoundingBoxes(2*j-1))], ...
                [min(segment_m(:,indexes.LaneCurvature)*2);min(segment_m(:,indexes.LaneCurvature)*2);max(segment_m(:,indexes.LaneCurvature)*2);max(segment_m(:,indexes.LaneCurvature)*2)], 'g', 'FaceAlpha', 0.2, 'LineStyle','none', "HandleVisibility","off");
        end
    end
    xlabel('X-UTM(m)'); ylabel("$\kappa(1/m)$");        
    set(gca,'FontSize', 14);
    xlim(relSegment);
    title("Road curvature");

    subplot(3,1,3);
    leftEdgeCoordinates = [X_abs-segment_m(:,indexes.LaneEdgePositionLeft).*sin(segment_m(:,indexes.theta_calc)) ...
        segment_m(:,indexes.Y_abs)+segment_m(:,indexes.LaneEdgePositionLeft).*cos(segment_m(:,indexes.theta_calc))];
    rightEdgeCoordinates = [X_abs-segment_m(:,indexes.LaneEdgePositionRight).*sin(segment_m(:,indexes.theta_calc)) ...
        segment_m(:,indexes.Y_abs)+segment_m(:,indexes.LaneEdgePositionRight).*cos(segment_m(:,indexes.theta_calc))];
    midLaneCoordinates = [X_abs-segment_m(:,indexes.c0).*sin(segment_m(:,indexes.theta_calc)) ...
        segment_m(:,indexes.Y_abs)+segment_m(:,indexes.c0).*cos(segment_m(:,indexes.theta_calc))];
    fill([leftEdgeCoordinates(:,1); flip(rightEdgeCoordinates(:,1))], [leftEdgeCoordinates(:,2); flip(rightEdgeCoordinates(:,2))], 'g');
    hold on;
    xlim(relSegment);
    for j=1:length(straightBoundingBoxes)/2
        if (j==1)
            fill([X_abs(straightBoundingBoxes(2*j-1)); X_abs(straightBoundingBoxes(2*j)); X_abs(straightBoundingBoxes(2*j)); X_abs(straightBoundingBoxes(2*j-1))], ...
                [min(get(gca,'ylim'));min(get(gca,'ylim'));max(get(gca,'ylim'));max(get(gca,'ylim'))], 'g', 'FaceAlpha', 0.2, 'DisplayName', 'Straight sections', 'LineStyle','none');
        else
            fill([X_abs(straightBoundingBoxes(2*j-1)); X_abs(straightBoundingBoxes(2*j)); X_abs(straightBoundingBoxes(2*j)); X_abs(straightBoundingBoxes(2*j-1))], ...
                [min(get(gca,'ylim'));min(get(gca,'ylim'));max(get(gca,'ylim'));max(get(gca,'ylim'))], 'g', 'FaceAlpha', 0.2, 'LineStyle','none', "HandleVisibility","off");
        end
    end
    plot(X_abs, segment_m(:,indexes.Y_abs), 'color', 'k', 'LineWidth', 2, 'DisplayName', 'Reference');
    plot(midLaneCoordinates(:,1), midLaneCoordinates(:,2), 'color', 'k', 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', 'Centerline');
    grid on;
    xlabel("X-UTM(m)"); ylabel("Y-UTM(m)");
    
    annotation('textarrow',[0.7 0.65],[0.22 0.22], 'String','Driving direction', 'FontSize', 14);
    ylims = get(gca,'ylim');
    set(gca,'FontSize', 14);
    title("Vehicle path");
    ylim(ylims);
end

