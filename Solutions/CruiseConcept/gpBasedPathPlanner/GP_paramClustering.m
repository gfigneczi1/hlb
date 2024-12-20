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

f = figure('Position',[100 100 1200 480]);
set(f,'defaulttextInterpreter','latex') ;
set(f, 'defaultAxesTickLabelInterpreter','latex');  
set(f, 'defaultLegendInterpreter','latex');
%% Additional functionality: check traffic dependent offset selection
for i=1:size(segments.segments,2)-2
    segment = segments.segments(i).segment;
    name = segments.segments(i).name;
    [~, segment_m, indexes] = prepareInputForPlanner(segment);
    subplot(1,3,1);
    [N,X]= hist(-segment_m(segment_m(:,indexes.OncomingTrafficType)==0, indexes.c0), 20);
    KPI(i,:,1) = [mean(-segment_m(segment_m(:,indexes.OncomingTrafficType)==0, indexes.c0)), std(-segment_m(segment_m(:,indexes.OncomingTrafficType)==0, indexes.c0))];
    plot(X,N/numel(find(segment_m(:,indexes.OncomingTrafficType)==0)),'DisplayName', name(1:5), 'Marker', '.');
    hold on;
    xline(mean(-segment_m(segment_m(:,indexes.OncomingTrafficType)==0, indexes.c0)), 'LineWidth', 1.5, 'HandleVisibility','off');
    grid on;
    ylim([0,1]);
    xlim([-1,1]);

    subplot(1,3,2);
    [N,X]= hist(-segment_m(segment_m(:,indexes.OncomingTrafficType)==1, indexes.c0), 20);
    KPI(i,:,2) = [mean(-segment_m(segment_m(:,indexes.OncomingTrafficType)==1, indexes.c0)), std(-segment_m(segment_m(:,indexes.OncomingTrafficType)==1, indexes.c0))];
    plot(X,N/numel(find(segment_m(:,indexes.OncomingTrafficType)==1)),'DisplayName', name(1:5), 'Marker', '.');
    hold on;
    xline(mean(-segment_m(segment_m(:,indexes.OncomingTrafficType)==1, indexes.c0)), 'LineWidth', 1.5, 'HandleVisibility','off');
    grid on;
    ylim([0,1]);
    xlim([-1,1]);

    subplot(1,3,3);
    [N,X]= hist(-segment_m(segment_m(:,indexes.OncomingTrafficType)>1, indexes.c0), 20);
    KPI(i,:,3) = [mean(-segment_m(segment_m(:,indexes.OncomingTrafficType)>1, indexes.c0)), std(-segment_m(segment_m(:,indexes.OncomingTrafficType)>1, indexes.c0))];
    plot(X,N/numel(find(segment_m(:,indexes.OncomingTrafficType)>1)),'DisplayName', name(1:5), 'Marker', '.');
    hold on;
    xline(mean(-segment_m(segment_m(:,indexes.OncomingTrafficType)>1, indexes.c0)), 'LineWidth', 1.5, 'HandleVisibility','off');
    grid on;
    ylim([0,1]);
    xlim([-1,1]);    
end
subplot(1,3,1); legend; title("$p(\delta| o_t=0)$");
subplot(1,3,2); legend; title("$p(\delta| o_t=1)$");
subplot(1,3,3); legend; title("$p(\delta| o_t>1)$");

f2 = figure('Position',[100 100 600 480]);
set(f2,'defaulttextInterpreter','latex') ;
set(f2, 'defaultAxesTickLabelInterpreter','latex');  
set(f2, 'defaultLegendInterpreter','latex');

confidencePoints = [1:1:numel(KPI(:,1,1)) flip(1:1:numel(KPI(:,1,1)))];
confidenceBounds = [KPI(:,1,1)+2*sqrt(KPI(:,2,1)); KPI(:,1,1)-2*sqrt(KPI(:,2,1))];
fill(confidencePoints', confidenceBounds, 'b', 'DisplayName', '95\% confidence', 'FaceAlpha',0.1);

hold on

confidencePoints = [1:1:numel(KPI(:,1,2)) flip(1:1:numel(KPI(:,1,2)))];
confidenceBounds = [KPI(:,1,2)+2*sqrt(KPI(:,2,2)); KPI(:,1,2)-2*sqrt(KPI(:,2,2))];
fill(confidencePoints', confidenceBounds, 'm', 'DisplayName', '95\% confidence', 'FaceAlpha',0.1);

confidencePoints = [1:1:numel(KPI(:,1,3)) flip(1:1:numel(KPI(:,1,3)))];
confidenceBounds = [KPI(:,1,3)+2*sqrt(KPI(:,2,3)); KPI(:,1,3)-2*sqrt(KPI(:,2,3))];
fill(confidencePoints', confidenceBounds, 'r', 'DisplayName', '95\% confidence', 'FaceAlpha',0.1);

plot(KPI(:,1,1), 'b', 'DisplayName', 'Free driving', 'LineWidth', 2);
plot(KPI(:,1,2), 'm', 'DisplayName', 'Small oncoming vehicle', 'LineWidth', 2)
plot(KPI(:,1,3), 'r', 'DisplayName', 'Truck oncoming vehicle', 'LineWidth', 2)
grid on
set(gca,"FontSize", 14);
legend('Location','best', 'FontSize',10);
xlabel('Driver ID'); ylabel('$\delta(m)$');


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
    %pathToParams = "C:\git\KDP\publications\GP\results\sparseGP_chosenParams\kpis";
    pathToParams = "C:\git\KDP\publications\GP\results\sparseGP_fixedAcceleration\12drivers\kpis";
    paramGP = dir(fullfile(pathToParams,"ETA_*"));
    for fileID = 1:length(paramGP)
        paramDriverID = str2num(paramGP(fileID).name(strfind(paramGP(fileID).name, 'driver_')+7:strfind(paramGP(fileID).name, '.')-1));
        paramData = load(fullfile(paramGP(fileID).folder, paramGP(fileID).name));
        N = length(paramData.ETA(1).hyp_opt.cov)-1; % -1: sigma_f is excluded
        for npID=1:10
            GP_params((npID-1)*N+1:npID*N, paramDriverID) = (paramData.ETA(npID).hyp_opt.cov(1:end-1))'; % columns are the parameters of one driver
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

