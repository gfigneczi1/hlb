% This function is the cover function to call the GP sub-function for a
% given segment_m. segment_ms is an array of struct containing the driving data
% of multiple drivers. Config is a metadata struct containing information
% about e.g., the plot folders
close all; clear;
load('C:\database\KDP_HLB_GP\Dr008_Dr023_input_GP.mat');
config.root = "./";

MAXIMUM_INPUT_SIZE = 15000; % before snipetting and normalization
MAXIMUM_LENGTH_OF_SNIPPET = 15000; % after cutting, snippets will be of this length
MAXIMUM_NUMBER_OF_SNIPPET = 1; % number of snippets after cutting the preprocessed data
GENERATE_KERNEL_CORRELATION_PLOTS = false;
KERNEL_TYPE = "{'covSEard'}"; %"{'covSum',{{'covProd', {'covLINard', 'covLINard'}}, 'covSEard'}}";
RATIO_OF_TRAIN_VS_TOTAL = 0.7;
GENERATE_OFFSET_TIME_PLOTS = false;
MAXIMUM_PREVIEW_DISTANCE = 150; % from within preview information is extracted, unit: m
OUTPUT_STEP_DISTANCE = 10; % in meters
NUMBER_OF_PREVIEW_INFORMATION = 2; % maximum number
OUTPUT_SHIFT = [10, 39, 136]; %linspace(15,MAXIMUM_PREVIEW_DISTANCE,10); %10:OUTPUT_STEP_DISTANCE:MAXIMUM_PREVIEW_DISTANCE; % preview distance is divided into sections where the output will be predicted
MERGE_DATA_TABLES = false;
DRIVER_ID_IF_NOT_MERGED = 1;
SIMPLIFICATION_RATIO = 1;
EPOCH_CALCULATION = true;
EPOCH_NUMBER = 10;
SPARSE_DATA = false;
GREEDY_REDUCTION = true;
FILTER_OUTPUT = true;
LDM_NP = [10, 39, 136];

usedInputs = ones(9,8);
I = -1*(eye(8)-1);
usedInputs(2:end,:) = I; % take out one input at a time

%usedInputs(1,2:3) = 0; % no O_t no fO_t

kernels = ["@covSEard", ...
    "{'covSum', {'covSEard', 'covLINard'}}",  ...
    "{'covProd', {'covLINard', 'covLINard'}}", ...
    "{'covSum',{{'covProd', {'covLINard', 'covLINard'}}, 'covSEard'}}", ...
    "{'covSum', {'covLINard', 'covSEard', {'covPPard',3}}}", ...
    "{'covPPard',3}", ...
    "{'covSum', {'covLINard', {'covPPard',3}}}", ...
    "{'covSum', {'covSEard', {'covPPard',3}}}"];
%kernels = ["{'covSum',{{'covProd', {'covLINard', 'covLINard'}}, 'covSEard'}}"];
%kernels = ["{'covSum',{{'covProd', {'covLINard', 'covLINard'}}, 'covSEard'}}"];

if (MERGE_DATA_TABLES)
    for i=1:length(segments.segments)
        % loop through segment_ms and concatenate each to the previous one
        if (i==1)
            segment=segments.segments(i).segment;
        else
            segment=[segment; segments.segments(i).segment];
        end
    end
    name = 'mergedDrivers';
else
    segment = segments.segments(DRIVER_ID_IF_NOT_MERGED).segment;
    name = segments.segments(DRIVER_ID_IF_NOT_MERGED).name;
end

% transforming the struct to array for lower calculation time
[~, segment_m, indexes] = prepareInputForPlanner(segment);

temp_folder_path = config.root;
plots_folder_name = 'plots';
set(0,'DefaultFigureVisible','off');

% PARAMETERS
shiftOnOutputSelection = OUTPUT_SHIFT;  %shift the offset in time (positive means shift forward)
p = RATIO_OF_TRAIN_VS_TOTAL; %percentage of evaluation data from entire dataset

for driverID = 1:size(segments.segments,2)
    DRIVER_ID_IF_NOT_MERGED = driverID;
    segment = segments.segments(DRIVER_ID_IF_NOT_MERGED).segment;
    name = segments.segments(DRIVER_ID_IF_NOT_MERGED).name;
    % transforming the struct to array for lower calculation time
    [~, segment_m, indexes] = prepareInputForPlanner(segment);
    dT = mean(diff(segment_m(:, indexes.Relative_time)));
    % generating the LDM output
    input = [segment_m(:, indexes.OncomingVehicleTimeToPass), ...
                segment_m(:, indexes.OncomingTrafficType), ...
                segment_m(:, indexes.FrontTrafficType), ...
                segment_m(:, indexes.VelocityX), ...
                movmean(segment_m(:, indexes.AccelerationX),20), ...
                movmean(segment_m(:, indexes.YawRate),20), ...
                movmean(segment_m(:, indexes.LaneCurvature), 20), ...
                movmean(segment_m(:, indexes.c3), 200)];

    [P_array, delta_0] = trainLDM(input, -movmean(segment_m(:,indexes.c0),180), LDM_NP);
    [estimationLDM] = resimulateTimeSequence (input, P_array, delta_0, LDM_NP);

    % now with different input sets, predict the error to the instinctive
    % path
for kernelID = 1:size(usedInputs,1)
    usedInput = usedInputs(kernelID,:);
    KERNEL_TYPE = kernels(1);
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

        input = input(:,usedInput==1);
        % resizing due to prior path
        input = input(1:size(estimationLDM{shiftID},1),:);
    
        variablesPool = ["$t_{pass}$", "$o_{type}$", "$fo_{type}$", "$v_x$", "$a_x$", "$\omega$", "$\kappa_{road}$", "$d\kappa$"];
        variables = variablesPool(usedInput==1);

        output = zeros(size(input,1),1);
        delta = -segment_m(:, indexes.c0);
        delta_filtered = movmean(delta,180);
        for shiftIDonOutput=1:size(input,1)
            if (shiftOnOutput(shiftIDonOutput) > size(input,1))
                break;
            else
                if (FILTER_OUTPUT)
                    output(shiftIDonOutput,1) = delta_filtered(shiftOnOutput(shiftIDonOutput));
                else
                    output(shiftIDonOutput,1) = delta(shiftOnOutput(shiftIDonOutput));
                end
            end
        end

        % resizing due to prior path
        output = output(1:size(estimationLDM{shiftID},1),:);

        % now subtracting the LDM path, in the corresponsing node point
        % (shiftID)
        outputLDM = estimationLDM{shiftID}; % LDM based offset
        output = output-outputLDM; % compensation

        % simplification with clustering
        if (SIMPLIFICATION_RATIO < 1)
            % long straight sections make the training difficult.
            % therefore, sections with low curvature and low yawrate values
            % are filtered out, only WHEN no object is present. 
            % SIMPLIFICATION RATIO gives how many percent of such points
            % are filtered out from the estimation data!
            relevantPoints = ((((abs(input(:,7)) > 0.2*max(abs(input(:,7)))) & ...
                (abs(input(:,8)) > 0.2*max(abs(input(:,8))))) ...
                | (input(:,2)>0)));
            irrelevantPoints = find(relevantPoints == 0);
            nIrrelevant = numel(find(irrelevantPoints));
            relevantPointsIndex = [find(relevantPoints==1); irrelevantPoints(randperm(nIrrelevant, floor(SIMPLIFICATION_RATIO*nIrrelevant)))];
    
            input = input(relevantPointsIndex,:);
            output = output(relevantPointsIndex,:);
            relevantPoints = all((abs(input(:,4:8)-mean(input(:,4:8)))<=2*std(input(:,4:8)))');
            input = input(relevantPoints,:);
            output = output(relevantPoints,:);
        end        

        % PLOTTING input - output plots for visual checks
        f = figure('units','normalized','outerposition',[0 0 1 1]);
        set(f,'defaulttextInterpreter','latex') ;
        set(f, 'defaultAxesTickLabelInterpreter','latex');  
        set(f, 'defaultLegendInterpreter','latex');
        subplot(size(input,2)+1,1,1);
        plot(output, 'color', 'k');
        grid on;
        
        for j=1:size(input,2)
            subplot(size(input,2)+1,1,(j+1));
            plot(input(:,j), 'color', 'b'); grid on;
            ylabel(variables(j));
        end            
        
        savefig(f, fullfile(temp_folder_path, plots_folder_name,...
                            strcat('GP_inputOutputs', name(1:end-4), '.fig')));
        saveas(f, fullfile(temp_folder_path, plots_folder_name,...
                        strcat('GP_inputOutputs_', name(1:end-4), '.png')));
        close(f);
    
    
        % SHUFFLE
        N = size(input,1);
        shuffledIndeces = randperm(N);
        input = input(shuffledIndeces,:);
        output = output(shuffledIndeces, :);
        outputLDM = outputLDM(shuffledIndeces,:);
        

        % LIMIT DATA IF NEEDED
        % this is done before norm and central
        input = input(1:min(size(input,1),MAXIMUM_INPUT_SIZE),:);
        output = output(1:min(size(input,1),MAXIMUM_INPUT_SIZE), :);
        outputLDM = outputLDM(1:min(size(input,1),MAXIMUM_INPUT_SIZE), :);

        % NORM AND CENTRAL
        [output, cout, s_out] = normalize(output);
        [input, cin, s_in] = normalize(input);
        % normalize the ldm output according to the output factors
        outputLDM = (outputLDM-cout)/s_out;
        
        % SNIPETTING
        for snippetID = 1:min(MAXIMUM_NUMBER_OF_SNIPPET,floor(N/MAXIMUM_LENGTH_OF_SNIPPET))
            snippets{snippetID,1} = input((snippetID-1)*MAXIMUM_LENGTH_OF_SNIPPET+1:snippetID*MAXIMUM_LENGTH_OF_SNIPPET,:);
            snippets{snippetID,2} = output((snippetID-1)*MAXIMUM_LENGTH_OF_SNIPPET+1:snippetID*MAXIMUM_LENGTH_OF_SNIPPET,:);
            snippets{snippetID,3} = outputLDM((snippetID-1)*MAXIMUM_LENGTH_OF_SNIPPET+1:snippetID*MAXIMUM_LENGTH_OF_SNIPPET,:);
        end

        % SIMULATION
        for i=1:min(MAXIMUM_NUMBER_OF_SNIPPET,length(snippets))
            Din = snippets{i,1};
            Dout = snippets{i,2};
            DoutLDM = snippets{i,3};
            M = size(Din,1);
            % EVALUATION / VALIDATION DATA SELECTION
            estimationData = 1:1:round(p*M);
            validationData = round(p*M)+1:1:M;
            input_estimation = Din(estimationData,:);
            output_estimation = Dout(estimationData,1);
            input_validation = Din(validationData,:);
            output_validation = Dout(validationData,:);
            outputLDM_estimation = DoutLDM(estimationData,1);
            outputLDM_validation = DoutLDM(validationData,1);

            % resampling estimation data
            k = 1;
            input_estimation = input_estimation(1:k:end,:);
            output_estimation = output_estimation(1:k:end,:);
            outputLDM_estimation = outputLDM_estimation(1:k:end,:);
    
            %% Define GP 
            meanfunc = [];       % Start with a zero mean prior
            eval(strcat('covfunc = ',KERNEL_TYPE));    % Squared Exponental covariance function
            % ID problem
            
            likfunc = @likGauss;    % Gaussian likelihood
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
                hyp = struct('mean', [], 'cov', ones(1,3*size(input_estimation,2)+1), 'lik', -1);
            elseif (KERNEL_TYPE == "{'covPPard',3}")
                hyp = struct('mean', [], 'cov', ones(1,size(input_estimation,2)+1), 'lik', -1);
            elseif (KERNEL_TYPE == "{'covSum', {'covLINard', 'covSEard', {'covPPard',3}}}") 
                hyp = struct('mean', [], 'cov', ones(1,3*size(input_estimation,2)+2), 'lik', -1);
            elseif (KERNEL_TYPE == "{'covSum', {'covSEiso', 'covLINiso'}}")
                hyp = struct('mean', [], 'cov', [0, 0, 0], 'lik', -1);
            elseif (KERNEL_TYPE == "{'covSum', {'covSEiso', 'covLINiso', {'covPERiso',{@covRQiso}}}}")
                hyp = struct('mean', [], 'cov', [0, 0, 0, 0, 0, 0, 0], 'lik', -1);
            elseif (KERNEL_TYPE == "{'covSum', {'covRQiso',{'covPERiso',{@covRQiso}}, {'covProd', {'covLINiso', 'covLINiso'}}}}")
                hyp = struct('mean', [], 'cov', [0, 0, 0, 0, 0, 0, 0, 0, 0], 'lik', -1);
            elseif (KERNEL_TYPE == "@covSEard" || KERNEL_TYPE == "{'covSEard'}")
                hyp = struct('mean', [], 'cov', ones(1,size(input_estimation,2)+1), 'lik', -1);
            elseif (KERNEL_TYPE == "@covLINard")
                hyp = struct('mean', [], 'cov', ones(1,size(input_estimation,2)), 'lik', -1);
            elseif (KERNEL_TYPE == "{'covSum', {'covSEard', 'covLINard'}}")
                hyp = struct('mean', [], 'cov', ones(1,2*size(input_estimation,2)+1), 'lik', -1);
            elseif (KERNEL_TYPE == "{'covProd', {'covLINard', 'covLINard'}}")
                hyp = struct('mean', [], 'cov', ones(1,2*size(input_estimation,2)), 'lik', -1);
            elseif (KERNEL_TYPE == "{'covSum',{{'covProd', {'covLINard', 'covLINard'}}, 'covSEard'}}")
                hyp = struct('mean', [], 'cov', ones(1,3*size(input_estimation,2)+1), 'lik', -1);
            elseif (KERNEL_TYPE == "{'covSum', {'covLINard', {'covPPard',3}}}")
                hyp = struct('mean', [], 'cov', ones(1,2*size(input_estimation,2)+1), 'lik', -1);
            elseif (KERNEL_TYPE == "{'covSum', {'covSEard', {'covPPard',3}}}")
                hyp = struct('mean', [], 'cov', ones(1,2*size(input_estimation,2)+2), 'lik', -1);
            else
                hyp = struct('mean', [], 'cov', [0, 0], 'lik', -1);
            end

            if (EPOCH_CALCULATION)
                % cut the estimation data to epochs, and loop through it for
                % training
                epochSize = floor(size(input_estimation,1)/EPOCH_NUMBER);
                for iter = 1:1
                    for epochID=1:EPOCH_NUMBER
                        epoch{epochID,1} = input_estimation((epochID-1)*epochSize+1:epochID*epochSize,:);
                        epoch{epochID,2} = output_estimation((epochID-1)*epochSize+1:epochID*epochSize,:);
                        epoch{epochID,3} = outputLDM_estimation((epochID-1)*epochSize+1:epochID*epochSize,:);
                        hypInduced = hyp;
                        if (SPARSE_DATA)
                            % generate initial value for induced inputs
                            %xu = linspace(-1,1,100);
                            %xu = xu'*ones(1,size(input_estimation,2));
                            %stepSize = round(size(input_estimation,1)/100);
                            %xu = epoch{epochID,1}(1:stepSize:end,:);
                            for inputID=1:size(input_estimation,2)
                                xu(:,inputID) = linspace(min(input_estimation(:,inputID)), max(input_estimation(:,inputID)), 1750)';
                            end
%                             k = kmeans(epoch{epochID,1}, 6100);
%                             for kId = 1:numel(unique(k))
%                                 xu(kId,:) = mean(epoch{epochID,1}(k==kId,:));
%                             end
                            covInduced = {'apxSparse',covfunc,xu};           % inducing points
                            inf = @infGaussLik;
                            infv  = @(varargin) inf(varargin{:},struct('s',0.0));           % VFE, opt.s = 0
                            hypInduced.xu = xu; % adding induced inputs
                            hypInducedOpt = minimize(hypInduced,@gp,-100, infv,meanfunc,covInduced,likfunc,epoch{epochID,1}, epoch{epochID,2}); 
                            uInduced{epochID} = hypInducedOpt.xu;
                            yInduced{epochID} = gp(hypInduced,@infGaussLik,meanfunc,covfunc,likfunc,epoch{epochID,1}, epoch{epochID,2}, hypInduced.xu);  
                        elseif (GREEDY_REDUCTION)
                            % greedy reduction: iteratively increase input
                            % data, then run optimization of
                            % hyperparameters
                            eps = 250;
                            if (epochID==1)
                                uGreedy = epoch{epochID,1};
                                yGreedy = epoch{epochID,2};
                                yGreedyLdm = epoch{epochID,3};
                            else
                                for sampleID = 2:size(epoch{epochID,1})
                                    eucledianNorm = norm(uGreedy-epoch{epochID,1}(sampleID,:));
                                    if (eucledianNorm > eps)
                                        uGreedy = [uGreedy; epoch{epochID,1}(sampleID,:)];
                                        yGreedy = [yGreedy; epoch{epochID,2}(sampleID,:)];
                                        yGreedyLdm = [yGreedyLdm; epoch{epochID,3}(sampleID,:)];
                                    end
                                end
                            end
                        else
                            hyp = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, epoch{epochID,1}, epoch{epochID,2}); % Optimize the marginal likelihood
                        end
                    end
                end
                if (SPARSE_DATA)
                    input_induced = []; output_induced = [];
                    for epochID=1:EPOCH_NUMBER
                        input_induced = [input_induced; uInduced{epochID}];
                        output_induced = [output_induced; yInduced{epochID}];
                    end
                    hyp_opt = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, input_induced, output_induced); % Optimize the marginal likelihood
                    % checking performance degradation
                    [estimationInduced, deviationInduced] = gp(hyp_opt, @infGaussLik, meanfunc, covfunc, likfunc, input_induced, output_induced, input_estimation);

                    %% KPI-s
                    % RMS calculation
                    RMSinduced = sqrt(sum((estimationInduced-output_estimation).^2)/size(estimationInduced,1));
                    % NRMS - normalization on scale
                    W = max(output_estimation) - min(output_estimation);
                    NRMS_W_induced = RMSinduced/W;
                    % NRMS - normalization on absolute maximum
                    M = max(abs(output_estimation));
                    NRMS_M_induced = RMSinduced/M;
                    % mean variance
                    RMS_DEV_induced = sqrt(sum(deviationInduced.^2)/size(deviationInduced,1));

                    fprintf("Induction accuracy:\n");
                    fprintf("RMS value is: %f\n", RMSinduced);
                    fprintf("NRMS value based on range: %f\n", NRMS_W_induced);
                    fprintf("NRMS value based on absolute maximum: %f\n", NRMS_M_induced);

                    % assigning for further analysis
                    input_estimation = input_induced;
                    output_estimation = output_induced;
                elseif (GREEDY_REDUCTION)
                    hyp_opt = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, uGreedy, yGreedy); 

                    % checking performance degradation
                    [estimationgreedy, deviationgreedy] = gp(hyp_opt, @infGaussLik, meanfunc, covfunc, likfunc, uGreedy, yGreedy, input_estimation);

                    %% KPI-s
                    % RMS calculation
                    RMSgreedy = sqrt(sum((estimationgreedy+outputLDM_estimation-(output_estimation+outputLDM_estimation)).^2)/size(estimationgreedy,1));
                    % NRMS - normalization on scale
                    W = max(output_estimation+outputLDM_estimation) - min(output_estimation+outputLDM_estimation);
                    NRMS_W_greedy = RMSgreedy/W;
                    % NRMS - normalization on absolute maximum
                    M = max(abs(output_estimation+outputLDM_estimation));
                    NRMS_M_greedy = RMSgreedy/M;
                    % mean variance
                    RMS_DEV_greedy = sqrt(sum(deviationgreedy.^2)/size(deviationgreedy,1));

                    fprintf("Induction accuracy:\n");
                    fprintf("RMS value is: %f\n", RMSgreedy);
                    fprintf("NRMS value based on range: %f\n", NRMS_W_greedy);
                    fprintf("NRMS value based on absolute maximum: %f\n", NRMS_M_greedy);

                    % assigning for further analysis
                    input_estimation = uGreedy;
                    output_estimation = yGreedy;
                    outputLDM_estimation = yGreedyLdm;
                else
                    hyp_opt = hyp;
                end
            else
                hyp_opt = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, input_estimation, output_estimation); % Optimize the marginal likelihood
            end
    
            if (RATIO_OF_TRAIN_VS_TOTAL < 1)
                % Evaluation of validation data
                [estimation, deviation] = gp(hyp_opt, @infGaussLik, meanfunc, covfunc, likfunc, input_estimation, output_estimation, input_validation); % extract the mean and covarance functions
                %% KPI-s
                % RMS calculation
                RMS = sqrt(sum((estimation+outputLDM_validation-(output_validation+outputLDM_validation)).^2)/size(estimation,1));
                % NRMS - normalization on scale
                W = max(output_validation+outputLDM_validation) - min(output_validation+outputLDM_validation);
                NRMS_W = RMS/W;
                % NRMS - normalization on absolute maximum
                M = max(abs(output_validation+outputLDM_validation));
                NRMS_M = RMS/M;
                % mean variance
                RMS_DEV = sqrt(sum(deviation.^2)/size(deviation,1));
                
                fprintf("EVALUTATION of snippet %d:\n", i);
                fprintf("RMS value is: %f\n", RMS);
                fprintf("NRMS value based on range: %f\n", NRMS_W);
                fprintf("NRMS value based on absolute maximum: %f\n", NRMS_M);
            else
                RMS=0; NRMS_W = 0; NRMS_M = 0; RMS_DEV = 0;
            end
    
            % Evaluation of estimation data
            [estimationEstimation, deviationEstimation] = gp(hyp_opt, @infGaussLik, meanfunc, covfunc, likfunc, input_estimation, output_estimation, input_estimation); % extract the mean and covarance functions
    
            %% KPI-s
            % RMS calculation
            RMSest = sqrt(sum((estimationEstimation+outputLDM_estimation-(output_estimation+outputLDM_estimation)).^2)/size(estimationEstimation,1));
            % NRMS - normalization on scale
            W = max(output_estimation+outputLDM_estimation) - min(output_estimation+outputLDM_estimation);
            NRMS_W_est = RMSest/W;
            % NRMS - normalization on absolute maximum
            M = max(abs(output_estimation+outputLDM_estimation));
            NRMS_M_est = RMSest/M;
            % mean variance
            RMS_DEV_est = sqrt(sum(deviationEstimation.^2)/size(deviationEstimation,1));
            
            fprintf("EVALUTATION of snippet %d:\n", i);
            fprintf("RMS_est value is: %f\n", RMSest);
            fprintf("NRMS_est value based on range: %f\n", NRMS_W_est);
            fprintf("NRMS_est value based on absolute maximum: %f\n", NRMS_M_est);
    
            f = figure('units','normalized','outerposition',[0 0 1 1]);
        
            if (RATIO_OF_TRAIN_VS_TOTAL < 1)
                % plot the estimation data
                subplot(2,1,1);
                confidenceBounds = [estimation+outputLDM_validation+2*sqrt(deviation); flip(estimation+outputLDM_validation-2*sqrt(deviation),1)];
                confidencePoints = (1:1:numel(estimation))';
            
                fill([confidencePoints; flip(confidencePoints)], confidenceBounds, 'y', 'DisplayName', '95% confidence');
                hold on;
                plot(estimation+outputLDM_validation,'LineWidth', 2, 'color', 'k', 'LineStyle','--', 'DisplayName', 'estimated by GP');    
                ylabel('offset');
                grid on;
                plot(output_validation+outputLDM_validation, 'LineWidth', 2, 'LineStyle', ':', 'color', 'b', 'DisplayName','validation data');
                legend;
                ylabel('offset(m)');
            end
    
            % plot the estimation data
            subplot(2,1,2);
            confidenceBounds = [estimationEstimation+outputLDM_estimation+2*sqrt(deviationEstimation); flip(estimationEstimation+outputLDM_estimation-2*sqrt(deviationEstimation),1)];
            confidencePoints = (1:1:numel(estimationEstimation))';
        
            fill([confidencePoints; flip(confidencePoints)], confidenceBounds, 'y', 'DisplayName', '95% confidence');
            hold on;
            plot(estimationEstimation+outputLDM_estimation,'LineWidth', 2, 'color', 'k', 'LineStyle','--', 'DisplayName', 'estimated by GP');    
            ylabel('offset');
            grid on;
            plot(output_estimation+outputLDM_estimation, 'LineWidth', 2, 'LineStyle', ':', 'color', 'b', 'DisplayName','estimation data');
            legend;
            ylabel('offset(m)');
            
            savefig(f, fullfile(temp_folder_path, plots_folder_name,...
                                strcat('Snippets_', num2str(i), '_', num2str(shiftID), '_driver', num2str(driverID), '.fig')));
            saveas(f, fullfile(temp_folder_path, plots_folder_name,...
                            strcat('Snippets_', num2str(i), '_', num2str(shiftID), '_driver', num2str(driverID), '.png')));
            close(f);

            if (GENERATE_KERNEL_CORRELATION_PLOTS && i==1 && shiftID==1)
            % Evaluation of entire range - Kernel shape check
            for j=1:size(input_estimation,2)
                input_range(:,j) = linspace(-1, 1, 30)'; 
            end
                    
            for k=1:size(input_range,2)
                f = figure('units','normalized','outerposition',[0 0 1 1]);
                set(f,'defaulttextInterpreter','latex') ;
                set(f, 'defaultAxesTickLabelInterpreter','latex');  
                set(f, 'defaultLegendInterpreter','latex');
                fprintf("Kernel generation k=%d/%d\n", k, size(input_range,2));
                for n=k+1:size(input_range,2)
                    for j=1:size(input_range,1)
                        input_estimation_range = input_range;
                        input_estimation_range(:,k) = input_range(j,k);
                        nonIndexed = 1:1:size(input_range,2);
                        input_estimation_range(:,nonIndexed~=k & nonIndexed~=n) = 0;
                        [estimationRange, deviationRange] = gp(hyp_opt, @infGaussLik, meanfunc, covfunc, likfunc, input_estimation, output_estimation, input_estimation_range); % extract the mean and covarance functions
    
                        estimationRangeArray(k,n, :,j) = estimationRange;
                        % the j^th variable in the k^th column is constant
                        % during this for cycle, which means one column of
                        % estimationRangeArray(k,n) reflect to the change in
                        % the 2nd variable, not the k^th variable (but the
                        % n^th).
                    end
                    subplot(2,4,n);
                    a(:,:) = estimationRangeArray(k,n,:,:);
                    surf(input_range(:,k), input_range(:,n), a);
                    hold on;
                    %plot3(input_estimation(:,k), input_estimation(:,n), output_estimation(:,1), 'kx');
                    xlabel(variables(k)); ylabel(variables(n)); zlabel("$\delta$");
                    set(gca,'FontSize', 14);
                end
                savefig(f, fullfile(temp_folder_path, plots_folder_name,...
                                strcat('Kernel_', num2str(k), '_', num2str(i), '_driver_', num2str(driverID), '_', name(1:end-4), '.fig')));
                saveas(f, fullfile(temp_folder_path, plots_folder_name,...
                                strcat('Kernel_', num2str(k), '_', num2str(i),'_driver_', num2str(driverID), '_', name(1:end-4), '.png')));
                close(f);
            end

            end
           
    
            if (GENERATE_OFFSET_TIME_PLOTS)
                [estimationEstimation, deviationEstimation] = gp(hyp_opt, @infGaussLik, meanfunc, covfunc, likfunc, input_estimation, output_estimation, input);
                f = figure();
                confidenceBounds = [estimationEstimation+2*sqrt(deviationEstimation); flip(estimationEstimation-2*sqrt(deviationEstimation),1)];
                confidencePoints = (1:1:numel(estimationEstimation))';
            
                fill([confidencePoints; flip(confidencePoints)], confidenceBounds, 'y', 'DisplayName', '95% confidence');
                hold on;
                plot(estimationEstimation,'LineWidth', 2, 'color', 'k', 'LineStyle','--', 'DisplayName', 'estimated by GP');    
                ylabel('offset');
                grid on;
                plot(output, 'LineWidth', 2, 'LineStyle', ':', 'color', 'b', 'DisplayName','validation data');
                legend;
                ylabel('offset(m)');
                xlabel('samples');
        
                savefig(f, fullfile(temp_folder_path, plots_folder_name,...
                                    strcat('TimePlot_', num2str(i), name(1:end-4), '.fig')));
                saveas(f, fullfile(temp_folder_path, plots_folder_name,...
                                strcat('TimePlot_', num2str(i), name(1:end-4), '.png')));
                close(f);
            end
        
            KPI{shiftID,i}= [RMS NRMS_W NRMS_M RMS_DEV RMSest NRMS_W_est NRMS_M_est RMS_DEV_est hyp_opt.cov];
            ETA(shiftID).hyp_opt = hyp_opt.cov;
            ETA(shiftID).input_estimation = input_estimation;
            ETA(shiftID).output_estimation = output_estimation;
            ETA(shiftID).output_estimationLDM = outputLDM_estimation;

            fprintf("time of shift id %d is %f\n", shiftID, toc);        
        end
    
        KPIsum(shiftID,:) = [mean(KPI{1}(:,2)) max(KPI{1}(:,2)) std(KPI{1}(:,3))];
    
    end % end of shift selection
    save( fullfile(temp_folder_path, plots_folder_name,strcat('KPI_input_', num2str(kernelID), '_driver_', num2str(driverID), '.mat')), 'KPI');
    save( fullfile(temp_folder_path, plots_folder_name,strcat('ETA_input_', num2str(kernelID), '_driver_', num2str(driverID), '.mat')), 'ETA');
    save( fullfile(temp_folder_path, plots_folder_name,strcat('KPIsum_input_', num2str(kernelID), '_driver_', num2str(driverID), '.mat')), 'KPIsum');
end
end

function curveTypes = calculateCurveType(segment_m_m, indexes, name)

global temp_folder_path plots_folder_name
    % return: curveTypes: 0 = unknown, 1=straight, 2=left, 3=right
    thd = 3.5e-4;
    curveTypes = zeros(size(segment_m_m,1),1);
    curvature = movmean(segment_m_m(:,indexes.c2)*2, 50);
    straightLine = abs(curvature) < thd;
    straightLine = morphologyOpen(straightLine, 100);
    straightLine = morphologyOpen(straightLine, 50);
    straightLine = morphologyClose(straightLine, 100);
    curveTypes(straightLine) = 1;
    
    leftCurve = curvature>=thd;
    leftCurve = morphologyClose(leftCurve,50);
    rightCurve = curvature<=-thd;
    rightCurve = morphologyClose(rightCurve,50);
    curveTransition = leftCurve&rightCurve;
    curveTypes(leftCurve&~curveTransition) = 2;
    curveTypes(rightCurve&~curveTransition) = 3;
    curveTypes(curveTransition) = 2.5;
end

function dataOut = morphologyOpen(dataIn, windowSize)
% this function is a morphology open filter with window size filter
dataOut = dataIn;
for i=1:length(dataIn)-windowSize
    if (~all(dataIn(i:i+windowSize)))
        dataOut(i:i+windowSize) = 0;
    end
end

end

function dataOut = morphologyClose(dataIn, windowSize)
% this function is a morphology close filter with window size filter
dataOut = dataIn;
for i=1:length(dataIn)-windowSize
    if (any(dataIn(i:i+windowSize)))
        dataOut(i:i+windowSize) = 1;
    end
end

end

function inputVariations = generateInputVariations(numberOfInputs)

x = 1:2^numberOfInputs-1;
inputVariations = dec2bin(x', numberOfInputs) == '1';

end

function [input,c,s] = normAndCentral(input)
for i=1:size(input,2)
     [input(:,i), c(i), s(i)] = normalize(input(:,i));
end
end

function [estimationLDM] = resimulateTimeSequence (input, P_array, delta_0, LDM_NP)
    % LDM resimulate
    dT = 0.05;
    ds = input(:,4)*dT;
    for i=1:size(input,1)
        S = cumsum(ds(i:end));
        np1 = i+find(S>=LDM_NP(1),1);
        np2 = i+find(S>=LDM_NP(2),1);
        np3 = i+find(S>=LDM_NP(3),1)-1;
        if (isempty(np3))
            break;
        end
        input_loc(i,:) = [mean(input(i:np1,7)) mean(input(np1:np2,7)) mean(input(np2:np3,7))];
    end
    leftCurves = mean(input_loc') >= 2.5e-4;
    rightCurves = mean(input_loc') <= -2.5e-4;
    inputLeft = input_loc; inputLeft(~leftCurves,:) = 0;
    inputRight = input_loc; inputRight(~rightCurves,:) = 0;

    for i=1:length(P_array)/2
        estimationLDM{i} = inputLeft/0.001*P_array{i}+inputRight/0.001*P_array{i+3}+delta_0;
    end
end


function [P_array, delta_0] = trainLDM(input, output, LDM_NP)
% LDM requires special input set, therefore a local input_loc is created 
% with necessary data
dT = 0.05;
ds = input(:,4)*dT;
delta_0 = mean(output(:,1));
for i=1:size(input,1)
    S = cumsum(ds(i:end));
    np1 = i+find(S>=LDM_NP(1),1);
    np2 = i+find(S>=LDM_NP(2),1);
    np3 = i+find(S>=LDM_NP(3),1)-1;
    if (isempty(np3))
        break;
    end
    input_loc(i,:) = [mean(input(i:np1,7)) mean(input(np1:np2,7)) mean(input(np2:np3,7))];
    output_loc(i,:) = output([np1, np2, np3],1)-delta_0;
end
P = functional_driverModelLearning(input_loc/0.001, output_loc ,8);
P = reshape(P(1:18),3,6);

P_array{1} = P(:,1);
P_array{2} = P(:,2);
P_array{3} = P(:,3);
P_array{4} = P(:,4);
P_array{5} = P(:,5);
P_array{6} = P(:,6);

end

