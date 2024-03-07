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
GENERATE_OFFSET_TIME_PLOTS = true;
MAXIMUM_PREVIEW_DISTANCE = 140; % from within preview information is extracted, unit: m
OUTPUT_STEP_DISTANCE = 10; % in meters
NUMBER_OF_PREVIEW_INFORMATION = 1; % maximum number
OUTPUT_SHIFT = 0; %10:OUTPUT_STEP_DISTANCE:MAXIMUM_PREVIEW_DISTANCE; % preview distance is divided into sections where the output will be predicted

segment = segments.segments(2).segment; % Driver 13
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
shiftOnOutputSelection = OUTPUT_SHIFT;  %shift the offset in time (positive means shift forward)
p = RATIO_OF_TRAIN_VS_TOTAL; %percentage of evaluation data from entire dataset

inputVariations = generateInputVariations(8);
inputVariations = ones(1,8);

dT = mean(diff(snippets{i}(:, indexes.q_T0)));

for numberOfPreviewInformation=1:NUMBER_OF_PREVIEW_INFORMATION
    previewSegmentLength = MAXIMUM_PREVIEW_DISTANCE/numberOfPreviewInformation; % unit: m
    %inputVariations = ones(1, 8+numberOfPreviewInformation); inputVariations(1) = 0;

for missingInputID=1:size(inputVariations,1)

for shiftID=1:numel(OUTPUT_SHIFT)
    for i=1:min(20,length(snippets))
        dx = snippets{i}(:, indexes.velocityX)*dT;

        % sweep a selection of shifts
        shiftOnOutput = [1:1:size(snippets{i}(:,indexes.q_T0),1)]'+floor(OUTPUT_SHIFT(shiftID)./dx);

        inputPreview = [];
        % input preview set is currently the mean curvature
        for previewInformationID = 1:numberOfPreviewInformation            
            if (numberOfPreviewInformation==1)
            else
                previewSegmentFrom = [1:1:size(snippets{i}(:,indexes.q_T0),1)]'+floor((previewSegmentLength*(previewInformationID-1))./dx);
                previewSegmentTo = [1:1:size(snippets{i}(:,indexes.q_T0),1)]'+floor((previewSegmentLength*previewInformationID)./dx);
                for dataPoint = 1:size(snippets{i}(:,indexes.q_T0),1)
                    if (previewSegmentTo(dataPoint) > size(snippets{i}(:,indexes.q_T0),1))
                        break;
                    else
                        inputPreview(dataPoint,previewInformationID) = mean(snippets{i}(previewSegmentFrom(dataPoint):previewSegmentTo(dataPoint), indexes.c2)*2);
                    end
                end
            end
        end

        %Npreview = size(inputPreview,1);
                
        input = [snippets{i}(:, indexes.c2)*2, ...
            snippets{i}(:, indexes.timeToPass), ...
            snippets{i}(:, indexes.oncomingTraffic), ...
            -movmean(snippets{i}(:, indexes.c1),20), ...
            snippets{i}(:, indexes.velocityX).*cos(-atan(movmean(snippets{i}(:, indexes.c1),20))), ...
            snippets{i}(:, indexes.velocityX).*sin(-atan(movmean(snippets{i}(:, indexes.c1),20))), ...
            movmean(snippets{i}(:, indexes.accelerationX),20), ...
            movmean(snippets{i}(:, indexes.yawRate),20)];
        
        input = normAndCentral(input);

        %input = [input(1:Npreview,:), inputPreview];
        
        input = input(:, find(inputVariations(missingInputID,:)==1));
    
        N = size(input,1);
        output = zeros(N,1);

%         for shiftIDonOutput=1:N
%             if (shiftOnOutput(shiftIDonOutput) > N)
%                 break;
%             else
%                 output(shiftIDonOutput,1) = -snippets{i}(shiftOnOutput(shiftIDonOutput), indexes.c0);
%             end
%         end
    
        % simple generation of output shift
        output(1:N-shiftOnOutput) = -snippets{i}(shiftOnOutput+1:end, indexes.c0); %-movmean(snippets{i}(shiftOnOutput+1:end, indexes.c0),180); %minus offset due to coordinate system transformation (vehicle to lane vs lane to vehicle)
    
        output = normAndCentral(output);
        
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
            input_range(:,j) = linspace(input_min(j), input_max(j), 100)'; 
        end
        
        variables = ["$\kappa$", "$t_{pass}$", "$obj_{type}$", "$\Theta$", "$v_x$", "$v_y$", "$a_x$", "$\omega$"];
        
        for k=1:size(input_range,2)
            f = figure('units','normalized','outerposition',[0 0 1 1]);
            set(f,'defaulttextInterpreter','latex') ;
            set(f, 'defaultAxesTickLabelInterpreter','latex');  
            set(f, 'defaultLegendInterpreter','latex');
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
                % 1st and 2nd dimensions based on the above comment
%                 scatter3(input_range(:,k), ...
%                	input_range(:,n), ...
%                 estimationRange, ...
%                 1, ...
%                 estimationRange);
                surf(input_range(:,k), input_range(:,n), a);
                xlabel(variables(k)); ylabel(variables(n));
                set(gca,'FontSize', 14);
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
    
        KPI{numberOfPreviewInformation, shiftID}(i,:) = [RMS NRMS_W NRMS_M hyp_opt.cov];
    
    end

    KPIsum(missingInputID, :) = [mean(KPI{1}(:,2)) max(KPI{1}(:,2)) std(KPI{1}(:,3))];
    
    f = figure();
    subplot(3,1,1);
    plot(KPI{numberOfPreviewInformation,shiftID}(:,1), 'Marker','x', 'LineWidth', 1.5, 'color', 'b');
    grid on;
    title(strcat('RMS value for', {' '}, num2str(shiftOnOutput), {' '}, 'shift on output'));
    ylabel('offset(m)'); xlabel('epoch');
    
    subplot(3,1,2);
    plot(KPI{numberOfPreviewInformation,shiftID}(:,2), 'Marker','x', 'LineWidth', 1.5, 'color', 'b');
    grid on;
    title(strcat('NRMS_W value for', {' '}, num2str(shiftOnOutput), {' '}, 'shift on output'));
    xlabel('epoch');
    
    subplot(3,1,3);
    plot(KPI{numberOfPreviewInformation,shiftID}(:,3), 'Marker','x', 'LineWidth', 1.5, 'color', 'b');
    grid on;
    title(strcat('NRMS_M value for', {' '}, num2str(shiftOnOutput), {' '}, 'shift on output'));
    xlabel('epoch');
    
    savefig(f, fullfile(temp_folder_path, plots_folder_name,...
                            strcat('KPI_', num2str(shiftID), name(1:end-4), '.fig')));
    saveas(f, fullfile(temp_folder_path, plots_folder_name,...
                    strcat('KPI_', num2str(shiftID), name(1:end-4), '.png')));
    close(f);

end % end of shift selection

%% SUMMARY PLOTS OF KPIs
f = figure();

for i=1:size(KPI,2)
    subplot(3,1,1);
    plot(shiftOnOutputSelection(i), mean(KPI{numberOfPreviewInformation,i}(:,1)), 'bo');
    hold on; grid on;
    title('Mean RMS for various selection')
    xlabel('shift on output'); ylabel('offset RMS(m)');

    subplot(3,1,2);
    plot(shiftOnOutputSelection(i), mean(KPI{numberOfPreviewInformation,i}(:,2)), 'bo');
    hold on; grid on;
    title('Mean NRMS_W for various selection')
    xlabel('shift on output'); ylabel('offset NRMS');

    subplot(3,1,3);
    plot(shiftOnOutputSelection(i), mean(KPI{numberOfPreviewInformation,i}(:,3)), 'bo');
    hold on; grid on;
    title('Mean NRMS_M for various selection')
    xlabel('shift on output'); ylabel('offset NRMS');
end

savefig(f, fullfile(temp_folder_path, plots_folder_name,...
                            strcat('KPISummary_', name(1:end-4), '_', num2str(numberOfPreviewInformation), '.fig')));
saveas(f, fullfile(temp_folder_path, plots_folder_name,...
                strcat('KPISummary_', name(1:end-4), '_', num2str(numberOfPreviewInformation), '.png')));
close(f);

save( fullfile(temp_folder_path, plots_folder_name,'KPI.mat'), 'KPI');

end
save( fullfile(temp_folder_path, plots_folder_name,'KPIsum.mat'), 'KPIsum');
end

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

function input = normAndCentral(input)
for i=1:size(input,2)
    input(:,i) = (input(:,i)-mean(input(:,i)))/max(abs(input(:,i)));
end
end

