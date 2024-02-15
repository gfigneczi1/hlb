% This function is the cover function to call the GP sub-function for a
% given segment. Segments is an array of struct containing the driving data
% of multiple drivers. Config is a metadata struct containing information
% about e.g., the plot folders
close all; clear;
load('C:\database\KDP_HLB_GP\Dr009_Dr013.mat');
config.root = "./";

MAXIMUM_LENGTH_OF_SNIPPET = 2500;
GENERATE_KERNEL_CORRELATION_PLOTS = false;
KERNEL_TYPE = "{'covSum', {'covRQiso',{'covPERiso',{@covRQiso}}, {'covProd', {'covLINiso', 'covLINiso'}}}}";
CUTTING_OPTION = "total";
RATIO_OF_TRAIN_VS_TOTAL = 0.7;
GENERATE_OFFSET_TIME_PLOTS = false;

segment = segments.segments(2).segment;
name = segments.segments(1).name;
% transforming the struct to array for lower calculation time
[~, segment_m, indexes] = prepareInputForPlanner(segment);

temp_folder_path = config.root;
plots_folder_name = 'plots';
set(0,'DefaultFigureVisible','off');

%% Input  and Output matrix (row: signals, column: time steps)
curveTypes = calculateCurveType(segment_m, indexes, name);
switch CUTTING_OPTION 
    case "snippeting"
        % Segmenting data to snippets with continuation
        j = 1;
        snippetStart = 1;
        for i=1:length(curveTypes)-1
            if (curveTypes(i) == 1 && curveTypes(i+1)~=1)
                % there was a step in the data, cut it
                snippets{j} = segment_m(snippetStart:i, :);
                j = j + 1;
            elseif (curveTypes(i) ~= 1 && curveTypes(i+1)==1)
                snippetStart = i+1;
            end
        end
    case "straights"
        straights = segment_m(curveTypes==1, :);
        if (numel(find(curveTypes==1)) > MAXIMUM_LENGTH_OF_SNIPPET)
            numberOfSnippets = floor(size(straights,1) / MAXIMUM_LENGTH_OF_SNIPPET);
            for i=1:numberOfSnippets
                snippets{i} = straights((i-1)*MAXIMUM_LENGTH_OF_SNIPPET+1:i*MAXIMUM_LENGTH_OF_SNIPPET,:);
            end
        else
            snippets{1} = segment_m(curveTypes==1, :);
        end
    case "total"
        if (size(segment_m,1) > MAXIMUM_LENGTH_OF_SNIPPET)
            numberOfSnippets = floor(size(segment_m,1) / MAXIMUM_LENGTH_OF_SNIPPET);
            for i=1:numberOfSnippets
                snippets{i} = segment_m((i-1)*MAXIMUM_LENGTH_OF_SNIPPET+1:i*MAXIMUM_LENGTH_OF_SNIPPET,:);
            end
        else
            snippets{1} = segment_m(:, :);
        end
    case "curves"
        % Segmenting data to snippets with continuation
        j = 1;
        snippetStart = 1;
        for i=1:length(curveTypes)-1
            if (curveTypes(i) > 1 && curveTypes(i+1)<=1)
                % there was a step in the data, cut it
                snippets{j} = segment_m(snippetStart:i, :);
                j = j + 1;
            elseif (curveTypes(i) <= 1 && curveTypes(i+1)>1)
                snippetStart = i+1;
            end
        end
end

%% loop through data sections
% PARAMETERS
shiftOnOutputSelection = 0;  %shift the offset in time (positive means shift forward)
p = RATIO_OF_TRAIN_VS_TOTAL; %percentage of evaluation data from entire dataset
inputVariations = generateInputVariations(8);

for missingInputID=211:size(inputVariations,1)

for shiftID=1:numel(shiftOnOutputSelection)
    % sweep a selection of shifts
    shiftOnOutput = shiftOnOutputSelection(shiftID);

    for i=1:length(snippets)
    
        ayRel = snippets{i}(:, indexes.accelerationY)-  snippets{i}(:,indexes.c2)*2.*snippets{i}(:,indexes.velocityX).^2;
                
        input = [snippets{i}(:, indexes.c2)*2, ...
            snippets{i}(:, indexes.timeToPass), ...
            snippets{i}(:, indexes.oncomingTraffic), ...
            -movmean(snippets{i}(:, indexes.c1),20), ...
            snippets{i}(:, indexes.velocityX).*cos(-atan(movmean(snippets{i}(:, indexes.c1),20))), ...
            snippets{i}(:, indexes.velocityX).*sin(-atan(movmean(snippets{i}(:, indexes.c1),20))), ...
            movmean(snippets{i}(:, indexes.accelerationX),20), ...
            movmean(snippets{i}(:, indexes.yawRate),20)];
        
        input = input(:, find(inputVariations(missingInputID,:)==1));
    
        N = size(input,1);
        output = zeros(N,1);
    
        output(1:N-shiftOnOutput) = -movmean(snippets{i}(shiftOnOutput+1:end, indexes.c0),180); %minus offset due to coordinate system transformation (vehicle to lane vs lane to vehicle)
    
            % PLOTTING input - output plots for visual checks
            f = figure('units','normalized','outerposition',[0 0 1 1]);
            subplot(size(input,2)+1,1,1);
            plot(output, 'color', 'k');
            grid on;
            
            for j=1:size(input,2)
                subplot(size(input,2)+1,1,(j+1));
                plot(input(:,j), 'color', 'b'); grid on;
            end            
            
            savefig(f, fullfile(temp_folder_path, plots_folder_name,...
                                strcat('GP_inputOutputs', name(1:end-4), '_', num2str(i), '.fig')));
            saveas(f, fullfile(temp_folder_path, plots_folder_name,...
                            strcat('GP_inputOutputs_', name(1:end-4), '_', num2str(i), '.png')));
            close(f);
    
        % EVALUATION / VALIDATION DATA SELECTION
        shuffledIndeces = randperm(N);
        estimationData = shuffledIndeces(1:round(p*N));
        validationData = shuffledIndeces(round(p*N)+1:end);
        input_estimation = input(estimationData,:);
        output_estimation = output(estimationData,1);
        input_validation = input(validationData,:);
        output_validation = output(validationData,:);
    
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
        elseif (KERNEL_TYPE == "{'covSum', {'covRQiso',{'covPERiso',{@covRQiso}}}}")
            hyp = struct('mean', [], 'cov', [0, 0, 0, 0, 0, 0, 0], 'lik', -1);
        elseif (KERNEL_TYPE == "{'covSum', {'covRQiso',{'covPERiso',{@covRQiso}}, 'covLINiso'}}")
            hyp = struct('mean', [], 'cov', [0, 0, 0, 0, 0, 0, 0, 0], 'lik', -1);
        elseif (KERNEL_TYPE == "{'covSum', {'covRQiso','covLINiso'}}")
            hyp = struct('mean', [], 'cov', [0, 0, 0, 0], 'lik', -1);
        elseif (KERNEL_TYPE == "{'covSum', {'covSEiso', 'covLINiso'}}")
            hyp = struct('mean', [], 'cov', [0, 0, 0], 'lik', -1);
        elseif (KERNEL_TYPE == "{'covSum', {'covSEiso', 'covLINiso', {'covPERiso',{@covRQiso}}}}")
            hyp = struct('mean', [], 'cov', [0, 0, 0, 0, 0, 0, 0], 'lik', -1);
        elseif (KERNEL_TYPE == "{'covSum', {'covRQiso',{'covPERiso',{@covRQiso}}, {'covProd', {'covLINiso', 'covLINiso'}}}}")
            hyp = struct('mean', [], 'cov', [0, 0, 0, 0, 0, 0, 0, 0, 0], 'lik', -1);
        else
            hyp = struct('mean', [], 'cov', [0, 0], 'lik', -1);
        end

        hyp_opt = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, input_estimation, output_estimation); % Optimize the marginal likelihood
       
        % Evaluation of validation data
        [estimation, deviation] = gp(hyp_opt, @infGaussLik, meanfunc, covfunc, likfunc, input_estimation, output_estimation, input_validation); % extract the mean and covarance functions
        
        % Evaluation of estimation data
        [estimationEstimation, deviationEstimation] = gp(hyp_opt, @infGaussLik, meanfunc, covfunc, likfunc, input_estimation, output_estimation, input_estimation); % extract the mean and covarance functions

        if (GENERATE_KERNEL_CORRELATION_PLOTS && i==1 && shiftID==1)
        % Evaluation of entire range - Kernel shape check
        input_min = min(input_estimation); input_max = max(input_estimation); 
        for j=1:size(input_estimation,2)
            input_range(:,j) = linspace(input_min(j), input_max(j), 20)'; 
        end

        for k=1:size(input_range,2)
            f = figure('units','normalized','outerposition',[0 0 1 1]);
            for n=k:size(input_range,2)
                for j=1:size(input_range,1)
                    input_estimation_range = input_range;
                    input_estimation_range(:,k) = input_estimation_range(j,k);
                    [estimationRange, deviationRange] = gp(hyp_opt, @infGaussLik, meanfunc, covfunc, likfunc, [input_estimation(:,k), input_estimation(:, n)], output_estimation, [input_estimation_range(:,k), input_estimation_range(:,n)]); % extract the mean and covarance functions
                    estimationRangeArray(k,n, :,j) = estimationRange;
                    % the j^th variable in the k^th column is constant
                    % during this for cycle, which means one column of
                    % estimationRangeArray(k,n) reflect to the change in
                    % the 2nd variable, not the k^th variable (but the
                    % n^th).
                end
                subplot(2,4,n);
                a(:,:) = estimationRangeArray(k,n,:,:);
                % 1st and 2nd dimensions based on the above comment
                surf(input_range(:,k), input_range(:,n), a);
                xlabel(strcat('var', num2str(k))); ylabel(strcat('var', num2str(n)));
            end
            savefig(f, fullfile(temp_folder_path, plots_folder_name,...
                            strcat('Kernel_', num2str(k), '_', num2str(i), '_', name(1:end-4), '.fig')));
            saveas(f, fullfile(temp_folder_path, plots_folder_name,...
                            strcat('Kernel_', num2str(k), '_', num2str(i),'_', name(1:end-4), '.png')));
            close(f);
        end
        end

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
        
        savefig(f, fullfile(temp_folder_path, plots_folder_name,...
                            strcat('Snippets_', num2str(i), '_', num2str(shiftID), name(1:end-4), '.fig')));
        saveas(f, fullfile(temp_folder_path, plots_folder_name,...
                        strcat('Snippets_', num2str(i), '_', num2str(shiftID), name(1:end-4), '.png')));
        close(f);
    
        %% KPI-s
        % RMS calculation
        RMS = sqrt(sum((estimation-output_validation).^2)/size(estimation,1));
        % NRMS - normalization on scale
        W = max(output_validation) - min(output_validation);
        NRMS_W = RMS/W;
        % NRMS - normalization on absolute maximum
        M = max(abs(output_validation));
        NRMS_M = RMS/M;
        fprintf("EVALUTATION of snippet %d:\n", i);
        fprintf("RMS value is: %f\n", RMS);
        fprintf("NRMS value based on range: %f\n", NRMS_W);
        fprintf("NRMS value based on absolute maximum: %f\n", NRMS_M);
    
        KPI{shiftID}(i,:) = [RMS NRMS_W NRMS_M hyp_opt.cov];
    
    end

    KPIsum(missingInputID, :) = [mean(KPI{1}(:,2)) max(KPI{1}(:,2)) std(KPI{1}(:,3))];
    
    f = figure();
    subplot(3,1,1);
    plot(KPI{shiftID}(:,1), 'Marker','x', 'LineWidth', 1.5, 'color', 'b');
    grid on;
    title(strcat('RMS value for', {' '}, num2str(shiftOnOutput), {' '}, 'shift on output'));
    ylabel('offset(m)'); xlabel('epoch');
    
    subplot(3,1,2);
    plot(KPI{shiftID}(:,2), 'Marker','x', 'LineWidth', 1.5, 'color', 'b');
    grid on;
    title(strcat('NRMS_W value for', {' '}, num2str(shiftOnOutput), {' '}, 'shift on output'));
    xlabel('epoch');
    
    subplot(3,1,3);
    plot(KPI{shiftID}(:,3), 'Marker','x', 'LineWidth', 1.5, 'color', 'b');
    grid on;
    title(strcat('NRMS_M value for', {' '}, num2str(shiftOnOutput), {' '}, 'shift on output'));
    xlabel('epoch');
    
    savefig(f, fullfile(temp_folder_path, plots_folder_name,...
                            strcat('KPI_', num2str(shiftID), name(1:end-4), '.fig')));
    saveas(f, fullfile(temp_folder_path, plots_folder_name,...
                    strcat('KPI_', num2str(shiftID), name(1:end-4), '.png')));
    close(f);

    if (GENERATE_OFFSET_TIME_PLOTS)
        [estimationEstimation, deviationEstimation] = gp(hyp_opt, @infGaussLik, meanfunc, covfunc, likfunc, input_estimation, output_estimation, input);
    end

end % end of shift selection

%% SUMMARY PLOTS OF KPIs
f = figure();

for i=1:length(KPI)
    subplot(3,1,1);
    plot(shiftOnOutputSelection(i), mean(KPI{i}(:,1)), 'bo');
    hold on; grid on;
    title('Mean RMS for various selection')
    xlabel('shift on output'); ylabel('offset RMS(m)');

    subplot(3,1,2);
    plot(shiftOnOutputSelection(i), mean(KPI{i}(:,2)), 'bo');
    hold on; grid on;
    title('Mean NRMS_W for various selection')
    xlabel('shift on output'); ylabel('offset NRMS');

    subplot(3,1,3);
    plot(shiftOnOutputSelection(i), mean(KPI{i}(:,3)), 'bo');
    hold on; grid on;
    title('Mean NRMS_M for various selection')
    xlabel('shift on output'); ylabel('offset NRMS');
end

savefig(f, fullfile(temp_folder_path, plots_folder_name,...
                            strcat('KPISummary_', name(1:end-4), '.fig')));
saveas(f, fullfile(temp_folder_path, plots_folder_name,...
                strcat('KPISummary_', name(1:end-4), '.png')));
close(f);

save( fullfile(temp_folder_path, plots_folder_name,'KPI.mat'), 'KPI');

end
save( fullfile(temp_folder_path, plots_folder_name,'KPIsum.mat'), 'KPIsum');

function curveTypes = calculateCurveType(segment_m, indexes, name)

global temp_folder_path plots_folder_name
    % return: curveTypes: 0 = unknown, 1=straight, 2=left, 3=right
    thd = 3.5e-4;
    curveTypes = zeros(size(segment_m,1),1);
    curvature = movmean(segment_m(:,indexes.c2)*2, 50);
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

    g = figure();

    plot(segment_m(:,indexes.X_abs), segment_m(:, indexes.Y_abs));
    grid on; hold on;
    plot(segment_m(curveTypes==1, indexes.X_abs), segment_m(curveTypes==1, indexes.Y_abs), 'ko');
    legend('original', 'straights');
    savefig(g, fullfile(temp_folder_path, plots_folder_name,...
                strcat('Map', name(1:end-4), '.fig')));
    saveas(g, fullfile(temp_folder_path, plots_folder_name,...
                strcat('Map_', name(1:end-4), '.png')));
    close(g);

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


