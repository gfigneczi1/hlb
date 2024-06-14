% This function is the cover function to call the GP sub-function for a
% given segment_m. segment_ms is an array of struct containing the driving data
% of multiple drivers. Config is a metadata struct containing information
% about e.g., the plot folders
close all; clear;
load('C:\database\KDP_HLB_GP\Dr008_Dr023_input_GP.mat');
config.root = "./";
READ_GP_PARAMS = true;

PARAMS.MAXIMUM_INPUT_SIZE = 15000; % before snipetting and normalization
PARAMS.KERNEL_TYPE = "{'covSum', {'covSEard', {'covPPard',3}}}";
PARAMS.RATIO_OF_TRAIN_VS_TOTAL = 0.7;
PARAMS.MAXIMUM_PREVIEW_DISTANCE = 150; % from within preview information is extracted, unit: m
PARAMS.OUTPUT_STEP_DISTANCE = 10; % in meters
PARAMS.NUMBER_OF_PREVIEW_INFORMATION = 2; % maximum number
PARAMS.OUTPUT_SHIFT = linspace(15,PARAMS.MAXIMUM_PREVIEW_DISTANCE,10); %10:OUTPUT_STEP_DISTANCE:MAXIMUM_PREVIEW_DISTANCE; % preview distance is divided into sections where the output will be predicted
PARAMS.EPOCH_CALCULATION = true;
PARAMS.EPOCH_NUMBER = 10;
PARAMS.GREEDY_REDUCTION = true;
PARAMS.LDM_NP = [10, 39, 136];

temp_folder_path = config.root;
plots_folder_name = 'plots';
set(0,'DefaultFigureVisible','off');

% PARAMETERS
shiftOnOutputSelection = PARAMS.OUTPUT_SHIFT;  %shift the offset in time (positive means shift forward)
p = PARAMS.RATIO_OF_TRAIN_VS_TOTAL; %percentage of evaluation data from entire dataset

for driverID = 1:size(segments.segments,2)
    DRIVER_ID_IF_NOT_MERGED = driverID;
    segment = segments.segments(DRIVER_ID_IF_NOT_MERGED).segment;
    name = segments.segments(DRIVER_ID_IF_NOT_MERGED).name;
    % transforming the struct to array for lower calculation time
    [~, segment_m, indexes] = prepareInputForPlanner(segment);
    [input, output, inputRaw, outputRaw, c_out, s_out, c_in, s_in] = prepareData(segment_m, indexes, PARAMS);
    
    M = size(input,1);
    % EVALUATION / VALIDATION DATA SELECTION
    estimationData = 1:1:round(p*M);
    validationData = round(p*M)+1:1:M;
    input_estimation = input(estimationData,:);
    output_estimation = output(estimationData,:);
    input_validation = input(validationData,:);
    output_validation = output(validationData,:);
    
    if (READ_GP_PARAMS)
        % read the learnt data from previous runs
        pathToParams = "C:\git\KDP\publications\GP\results\gp_model";
        paramGP = dir(fullfile(pathToParams,"ETA_input*"));
        for fileID = 1:length(paramGP)
            paramDriverID = str2num(paramGP(fileID).name(strfind(paramGP(fileID).name, 'driver_')+7:strfind(paramGP(fileID).name, '.')-1));
            if (paramDriverID == driverID)
                paramData = load(fullfile(paramGP(fileID).folder, paramGP(fileID).name));
                for npID=1:10
                    GP_params.hyp_opt_array{npID} = paramData.ETA(npID).hyp_opt;
                    GP_params.output_estimation{npID} = paramData.ETA(npID).output_estimation;
                    GP_params.input_estimation{npID} = paramData.ETA(npID).input_estimation;
                end
                break;
            else
                paramData = [];
            end
        end
    else
        % Train GPs
        hyp_opt_array = trainGPs(input_estimation, output_estimation, PARAMS);
        GP_params.hyp_opt_array = hyp_opt_array;
        for i=1:10
            GP_params.input_estimation{i} = input_estimation;
            GP_params.output_estimation{i} = output_estimation(:,i);
        end
    end

    % Train LDM
    [LDM_params.P_array, LDM_params.delta_0] = trainLDM(inputRaw, outputRaw, PARAMS);
        
    % Train LRM
    P_array = trainLRM(input_estimation, output_estimation, PARAMS);
    LRM_params.P_array = P_array;
    
    % Generate output
    [inputResim, c, s] = normAndCentral(inputRaw);
    [estimationGP, deviationGP, estimationLRM, estimationLDM] = resimulateTimeSequence (inputResim, GP_params, LRM_params, LDM_params, PARAMS, c,s);

    for i=1:1
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
            straightBoundingBoxes(2:2:end) = find(diff(straightSections) < 0);
            straightBoundingBoxes(1:2:end) = [find(diff(straightSections) > 0)];
        else
            straightBoundingBoxes(1:2:end) = find(diff(straightSections) > 0);
            straightBoundingBoxes(2:2:end) = find(diff(straightSections) < 0);
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
        plot(X_abs, estimationLRM{i}*s_out(i)+c_out(i), 'color', 'b',  'LineWidth', 2, 'DisplayName', 'LRM');
        plot(X_abs(1:length(estimationLDM{i})), estimationLDM{i}, 'color', 'g',  'LineWidth', 2, 'DisplayName', 'LDM');
        plot(X_abs, outputRaw(:,i), 'color', 'k', 'LineStyle', '--', 'LineWidth', 2, 'DisplayName', 'Reference');
        xlabel('X-UTM(m)'); ylabel("$\delta(m)$");
        yline(0, "HandleVisibility","off", 'Alpha',0.3, 'Color','k', 'LineWidth',3);

        legend("Location", "best", "Orientation","horizontal");
        set(gca,'FontSize', 14);
        xlim(relSegment);
        ylim([-1,1]);
        title(strcat("Estimation of $\delta_1$ - Driver", {' '}, num2str(driverID)));

        subplot(3,1,2);
        plot(X_abs, segment_m(:,indexes.LaneCurvature)*2, 'LineWidth', 2, 'color', 'k');
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

        savefig(f, fullfile(temp_folder_path, plots_folder_name,...
                            strcat('Resimulate_plots_driver_', num2str(driverID), '.fig')));
        saveas(f, fullfile(temp_folder_path, plots_folder_name,...
                        strcat('Resimulate_plots_driver_', num2str(driverID), '.png')));
        close(f);
    end
    
end

%% SUPPORT FUNCTIONS
function [estimationGP, deviationGP, estimationLRM, estimationLDM] = resimulateTimeSequence (input_validation, GP_params, LRM_params, LDM_params, PARAMS, c,s)
    % GP resimulate
    meanfunc = [];       % Start with a zero mean prior
    eval(strcat('covfunc = ',PARAMS.KERNEL_TYPE));    % Squared Exponental covariance function
    % ID problem            
    likfunc = @likGauss;    % Gaussian likelihood
    for i = 1:length(GP_params.hyp_opt_array)
        hyp = struct('mean', [], 'cov', 0, 'lik', -1);
        hyp.cov = GP_params.hyp_opt_array{i};
        [estimationGP{i}, deviationGP{i}] = gp(hyp, @infGaussLik, meanfunc, covfunc, likfunc, GP_params.input_estimation{i}, GP_params.output_estimation{i}, input_validation); % extract the mean and covarance functions
    end
    
    % LRM resimulate
    for i=1:length(LRM_params.P_array)
        estimationLRM{i} = input_validation*LRM_params.P_array{i}';
    end

    % LDM resimulate
    % denormalize input data
    inputRaw = input_validation*diag(s)+c;
    dT = 0.05;
    ds = inputRaw(:,2)*dT;
    for i=1:size(inputRaw,1)
        S = cumsum(ds(i:end));
        np1 = i+find(S>=PARAMS.LDM_NP(1),1);
        np2 = i+find(S>=PARAMS.LDM_NP(2),1);
        np3 = i+find(S>=PARAMS.LDM_NP(3),1)-1;
        if (isempty(np3))
            break;
        end
        input_loc(i,:) = [mean(inputRaw(i:np1,5)) mean(inputRaw(np1:np2,5)) mean(inputRaw(np2:np3,5))];
    end
    for i=1:length(LDM_params.P_array)
        estimationLDM{i} = input_loc*LDM_params.P_array{i}'+LDM_params.delta_0;
    end
end

function [input, output, inputRaw, outputRaw, c_out, s_out, c_in, s_in] = prepareData(segment_m, indexes, PARAMS)
     % input array
     input = [segment_m(:, indexes.OncomingVehicleTimeToPass), ...
        segment_m(:, indexes.VelocityX), ...
        movmean(segment_m(:, indexes.AccelerationX),20), ...
        movmean(segment_m(:, indexes.YawRate),20), ...
        movmean(segment_m(:, indexes.LaneCurvature), 20), ...
        movmean(segment_m(:, indexes.c3), 200)];
    
    % output array
    output = zeros(size(input,1),numel(PARAMS.OUTPUT_SHIFT));
    for shiftID=1:numel(PARAMS.OUTPUT_SHIFT)
        dT = mean(diff(segment_m(:, indexes.Relative_time)));
        dx = segment_m(:, indexes.VelocityX)*dT;        
        shiftOnOutput = [1:1:size(segment_m(:,indexes.Relative_time),1)]'+floor(PARAMS.OUTPUT_SHIFT(shiftID)./dx);
        for shiftIDonOutput=1:size(input,1)
            if (shiftOnOutput(shiftIDonOutput) > size(input,1))
                break;
            else
                output(shiftIDonOutput,shiftID) = -segment_m(shiftOnOutput(shiftIDonOutput), indexes.c0);
            end
        end
    end
    
    inputRaw = input;
    outputRaw = output;
    
    % SHUFFLE
    N = size(input,1);
    shuffledIndeces = randperm(N);
    input = input(shuffledIndeces,:);
    output = output(shuffledIndeces, :);

    % LIMIT DATA IF NEEDED
    % this is done before norm and central
    input = input(1:min(size(input,1),PARAMS.MAXIMUM_INPUT_SIZE),:);
    output = output(1:min(size(input,1),PARAMS.MAXIMUM_INPUT_SIZE), :);

    % NORM AND CENTRAL
    [output, c_out, s_out] = normAndCentral(output);
    [input, c_in, s_in] = normAndCentral(input);
end

function [dataOut, c,s]= normAndCentral(dataIn)
    for i=1:size(dataIn,2)
        c(i) = mean(dataIn(:,i));
        s(i) = std(dataIn(:,i));
        dataOut(:,i) = (dataIn(:,i)-mean(dataIn(:,i)))/std(dataIn(:,i));
    end
end

function P_array = trainLRM(input, output, PARAMS)
    input_estimation = input;
    for shiftID=1:numel(PARAMS.OUTPUT_SHIFT)
        output_estimation = output(:,shiftID);
        P_array{shiftID} = output_estimation'*input_estimation*inv(input_estimation'*input_estimation);
    end
end

function [P_array, delta_0] = trainLDM(input, output, PARAMS)
% LDM requires special input set, therefore a local input_loc is created 
% with necessary data
dT = 0.05;
ds = input(:,2)*dT;
delta_0 = mean(output(:,1));
for i=1:size(input,1)
    S = cumsum(ds(i:end));
    np1 = i+find(S>=PARAMS.LDM_NP(1),1);
    np2 = i+find(S>=PARAMS.LDM_NP(2),1);
    np3 = i+find(S>=PARAMS.LDM_NP(3),1)-1;
    if (isempty(np3))
        break;
    end
    try
        input_loc(i,:) = [mean(input(i:np1,5)) mean(input(np1:np2,5)) mean(input(np2:np3,5))];
    catch
        i;
    end
    output_loc(i,:) = output([np1, np2, np3],1)-delta_0;
end
P_array{1} = output_loc(:,1)'*input_loc*inv(input_loc'*input_loc);
P_array{2} = output_loc(:,2)'*input_loc*inv(input_loc'*input_loc);
P_array{3} = output_loc(:,3)'*input_loc*inv(input_loc'*input_loc);

end


function hyp_opt_array = trainGPs(input, output, PARAMS)
    
    for shiftID=1:numel(PARAMS.OUTPUT_SHIFT)
        output_estimation = output(:,shiftID);
        input_estimation = input;
    
        %% Define GP 
        meanfunc = [];       % Start with a zero mean prior
        eval(strcat('covfunc = ',PARAMS.KERNEL_TYPE));    % Squared Exponental covariance function
        % ID problem
            
        likfunc = @likGauss;    % Gaussian likelihood
        if (PARAMS.KERNEL_TYPE == "@covLINiso")
            hyp = struct('mean', [], 'cov', 0, 'lik', -1); % initalize hyper-parameter structure associated with the mean, covariance an likelihood functions
        elseif (PARAMS.KERNEL_TYPE == "@covRQiso")
            hyp = struct('mean', [], 'cov', [0,0,0], 'lik', -1);
        elseif (PARAMS.KERNEL_TYPE=="{'covPERiso',{@covRQiso}}")
            hyp = struct('mean', [], 'cov', [0,0,0,0], 'lik', -1);
        elseif (PARAMS.KERNEL_TYPE == "{'covSum', {'covRQiso',{'covPERiso',{@covRQiso}}}}")
            hyp = struct('mean', [], 'cov', [0, 0, 0, 0, 0, 0, 0], 'lik', -1);
        elseif (PARAMS.KERNEL_TYPE == "{'covSum', {'covRQiso',{'covPERiso',{@covRQiso}}, 'covLINiso'}}")
            hyp = struct('mean', [], 'cov', [0, 0, 0, 0, 0, 0, 0, 0], 'lik', -1);
        elseif (PARAMS.KERNEL_TYPE == "{'covSum', {'covRQiso','covLINiso'}}")
            hyp = struct('mean', [], 'cov', [0, 0, 0, 0], 'lik', -1);
        elseif (PARAMS.KERNEL_TYPE == "{'covSum', {'covSEiso', 'covLINiso'}}")
            hyp = struct('mean', [], 'cov', [0, 0, 0], 'lik', -1);
        elseif (PARAMS.KERNEL_TYPE == "{'covSum', {'covSEiso', 'covLINiso', {'covPERiso',{@covRQiso}}}}")
            hyp = struct('mean', [], 'cov', [0, 0, 0, 0, 0, 0, 0], 'lik', -1);
        elseif (PARAMS.KERNEL_TYPE == "{'covSum', {'covRQiso',{'covPERiso',{@covRQiso}}, {'covProd', {'covLINiso', 'covLINiso'}}}}")
            hyp = struct('mean', [], 'cov', [0, 0, 0, 0, 0, 0, 0, 0, 0], 'lik', -1);
        elseif (PARAMS.KERNEL_TYPE == "@covSEard" || PARAMS.KERNEL_TYPE == "{'covSEard'}")
            hyp = struct('mean', [], 'cov', ones(1,size(input_estimation,2)+1), 'lik', -1);
        elseif (PARAMS.KERNEL_TYPE == "@covLINard")
            hyp = struct('mean', [], 'cov', ones(1,size(input_estimation,2)), 'lik', -1);
        elseif (PARAMS.KERNEL_TYPE == "{'covSum', {'covSEard', 'covLINard'}}")
            hyp = struct('mean', [], 'cov', ones(1,2*size(input_estimation,2)+1), 'lik', -1);
        elseif (PARAMS.KERNEL_TYPE == "{'covProd', {'covLINard', 'covLINard'}}")
            hyp = struct('mean', [], 'cov', ones(1,2*size(input_estimation,2)), 'lik', -1);
        elseif (PARAMS.KERNEL_TYPE == "{'covSum',{{'covProd', {'covLINard', 'covLINard'}}, 'covSEard'}}")
            hyp = struct('mean', [], 'cov', ones(1,3*size(input_estimation,2)+1), 'lik', -1);
        elseif (PARAMS.KERNEL_TYPE == "{'covSum', {'covSEard', {'covPPard',3}}}")
            hyp = struct('mean', [], 'cov', ones(1,2*size(input_estimation,2)+2), 'lik', -1);
        else
            hyp = struct('mean', [], 'cov', [0, 0], 'lik', -1);
        end

        % cut the estimation data to epochs, and loop through it for
        % training
        epochSize = floor(size(input_estimation,1)/PARAMS.EPOCH_NUMBER);
        for epochID=1:PARAMS.EPOCH_NUMBER
            epoch{epochID,1} = input_estimation((epochID-1)*epochSize+1:epochID*epochSize,:);
            epoch{epochID,2} = output_estimation((epochID-1)*epochSize+1:epochID*epochSize,:);
            hypInduced = hyp;
            if (PARAMS.GREEDY_REDUCTION)
                % greedy reduction: iteratively increase input
                % data, then run optimization of
                % hyperparameters
                eps = 175;
                if (epochID==1)
                    uGreedy = epoch{epochID,1};
                    yGreedy = epoch{epochID,2};
                else
                    for sampleID = 2:size(epoch{epochID,1})
                        eucledianNorm = norm(uGreedy-epoch{epochID,1}(sampleID,:));
                        if (eucledianNorm > eps)
                            uGreedy = [uGreedy; epoch{epochID,1}(sampleID,:)];
                            yGreedy = [yGreedy; epoch{epochID,2}(sampleID,:)];
                        end
                    end
                end
            else
                hyp = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, epoch{epochID,1}, epoch{epochID,2}); % Optimize the marginal likelihood
            end
        end

        if (PARAMS.GREEDY_REDUCTION)
            hyp_opt = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, uGreedy, yGreedy); 

            % assigning for further analysis
            input_estimation = uGreedy;
            output_estimation = yGreedy;
        else
            hyp_opt = hyp;
        end
        hyp_opt_array{shiftID} = hyp_opt;
        
        clear uGreedy yGreedy epoch
    end % end of shift selection
end

