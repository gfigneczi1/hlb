close all; clear;
load('C:\database\KDP_HLB_GP\Driver_with_traffic.mat');
config.root = "./";

MAXIMUM_LENGTH_OF_SNIPPET = 25000;
CUTTING_OPTION = "total";
RATIO_OF_TRAIN_VS_TOTAL = 0.7;
GENERATE_OFFSET_TIME_PLOTS = true;
MAXIMUM_PREVIEW_DISTANCE = 150; % from within preview information is extracted, unit: m
OUTPUT_STEP_DISTANCE = 10; % in meters
NUMBER_OF_PREVIEW_INFORMATION = 10; % maximum number
OUTPUT_SHIFT = linspace(0,MAXIMUM_PREVIEW_DISTANCE,10); %10:OUTPUT_STEP_DISTANCE:MAXIMUM_PREVIEW_DISTANCE; % preview distance is divided into sections where the output will be predicted
MERGE_DATA_TABLES = false;
DRIVER_ID_IF_NOT_MERGED = 5;
SIMPLIFICATION_RATIO = 0;

if (MERGE_DATA_TABLES)
    for i=1:length(segments.segments)
        % loop through segments and concatenate each to the previous one
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

for i=1:min(20,length(snippets))
    dT = mean(diff(snippets{i}(:, indexes.q_T0)));
    dx = snippets{i}(:, indexes.velocityX)*dT;
    
    input = [snippets{i}(:, indexes.timeToPass), ...
        snippets{i}(:, indexes.oncomingTraffic), ...
        -movmean(atan(snippets{i}(:, indexes.c1)),20), ...
        snippets{i}(:, indexes.velocityX).*cos(-atan(movmean(snippets{i}(:, indexes.c1),20))), ...
        snippets{i}(:, indexes.velocityX).*sin(-atan(movmean(snippets{i}(:, indexes.c1),20))), ...
        movmean(snippets{i}(:, indexes.accelerationX),20), ...
        movmean(snippets{i}(:, indexes.yawRate),20), ...
        movmean(snippets{i}(:, indexes.c2)*2, 20), ...
        movmean(snippets{i}(:, indexes.c3), 200)];

    input = normAndCentral(input);
    
    output = zeros(size(input,1),numel(OUTPUT_SHIFT));
    for shiftID=1:numel(OUTPUT_SHIFT)
        % sweep a selection of shifts
        shiftOnOutput = [1:1:size(snippets{i}(:,indexes.q_T0),1)]'+floor(OUTPUT_SHIFT(shiftID)./dx);
        for shiftIDonOutput=1:size(input,1)
            if (shiftOnOutput(shiftIDonOutput) > size(input,1))
                break;
            else
                output(shiftIDonOutput,shiftID) = -snippets{i}(shiftOnOutput(shiftIDonOutput), indexes.c0);
            end
        end
    end
    
    [output, mu, alfa] = normAndCentral(output);
    
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
                            strcat('LDM_inputOutputs', name(1:end-4), '_', num2str(i), '.fig')));
        saveas(f, fullfile(temp_folder_path, plots_folder_name,...
                        strcat('LDM_inputOutputs_', name(1:end-4), '_', num2str(i), '.png')));
        close(f);

    N = size(input,1);

    % EVALUATION / VALIDATION DATA SELECTION
    shuffledIndeces = randperm(N);
    estimationData = shuffledIndeces(1:round(p*N));
    validationData = shuffledIndeces(round(p*N)+1:end);
    input_estimation = input(estimationData,:);
    output_estimation = output(estimationData,:);
    input_validation = input(validationData,:);
    output_validation = output(validationData,:);

    % simplification with clustering
    if (SIMPLIFICATION_RATIO < 1)
        % long straight sections make the training difficult.
        % therefore, sections with low curvature and low yawrate values
        % are filtered out, only WHEN no object is present. 
        % SIMPLIFICATION RATIO gives how many percent of such points
        % are filtered out from the estimation data!
        relevantPoints = (((abs(input_estimation(:,7)) > 0.1) & ...
            (abs(input_estimation(:,8)) > 0.1)) ...
            | (input_estimation(:,2)>0));
        irrelevantPoints = find(relevantPoints == 0);
        nIrrelevant = numel(find(irrelevantPoints));
        relevantPointsIndex = [find(relevantPoints==1); irrelevantPoints(randperm(nIrrelevant, floor(SIMPLIFICATION_RATIO*nIrrelevant)))];
        input_estimation = input_estimation(relevantPointsIndex,:);
        output_estimation = output_estimation(relevantPointsIndex,:);
    end
    
    %% LDM regression
    P = output_estimation'*input_estimation*inv(input_estimation'*input_estimation);

    if (RATIO_OF_TRAIN_VS_TOTAL < 1)
        for shiftID = 1:size(output_estimation,2)
            % Evaluation of validation data
            estimation = input_validation*P';
            %% KPI-s
            % RMS calculation
            RMS(shiftID) = sqrt(sum((estimation(:,shiftID)-output_validation(:,shiftID)).^2)/size(estimation,1));
            % NRMS - normalization on scale
            W = max(output_validation(:,shiftID)) - min(output_validation(:,shiftID));
            NRMS_W(shiftID) = RMS(shiftID)/W;
            % NRMS - normalization on absolute maximum
            M = max(abs(output_validation(:,shiftID)));
            NRMS_M(shiftID) = RMS(shiftID)/M;
        end
        KPI(i,:) = [RMS NRMS_W NRMS_M];

        fprintf("EVALUTATION of snippet %d:\n", i);
        fprintf("RMS value is: %f\n", RMS);
        fprintf("NRMS value based on range: %f\n", NRMS_W);
        fprintf("NRMS value based on absolute maximum: %f\n", NRMS_M);
    end

    % Evaluation of estimation data
    estimationEstimation = input_estimation*P';
    
    %% KPI-s
    % RMS calculation
    for shiftID=1:size(output_estimation,2)
        RMSest(shiftID) = sqrt(sum((estimationEstimation(:, shiftID)-output_estimation(:,shiftID)).^2)/size(estimationEstimation,1));
        % NRMS - normalization on scale
        W = max(output_estimation(:,shiftID)) - min(output_estimation(:, shiftID));
        NRMS_W_est(shiftID) = RMSest(shiftID)/W;
        % NRMS - normalization on absolute maximum
        M = max(abs(output_estimation(:, shiftID)));
        NRMS_M_est(shiftID) = RMSest(shiftID)/M;
    end

    fprintf("EVALUTATION of snippet %d:\n", i);
    fprintf("RMS_est value is: %f\n", RMSest);
    fprintf("NRMS_est value based on range: %f\n", NRMS_W_est);
    fprintf("NRMS_est value based on absolute maximum: %f\n", NRMS_M_est);

    f = figure('units','normalized','outerposition',[0 0 1 1]);

    if (RATIO_OF_TRAIN_VS_TOTAL < 1)
        % plot the estimation data
        subplot(2,1,1);
        hold on;
        plot(estimation(:,1),'LineWidth', 2, 'color', 'k', 'LineStyle','--', 'DisplayName', 'estimated by GP');    
        ylabel('offset');
        grid on;
        plot(output_validation(:,1), 'LineWidth', 2, 'LineStyle', ':', 'color', 'b', 'DisplayName','validation data');
        legend;
        ylabel('offset(m)');
    end

    % plot the estimation data
    subplot(2,1,2);
    hold on;
    plot(estimationEstimation(:,1),'LineWidth', 2, 'color', 'k', 'LineStyle','--', 'DisplayName', 'estimated by LDM');    
    ylabel('offset');
    grid on;
    plot(output_estimation(:,1), 'LineWidth', 2, 'LineStyle', ':', 'color', 'b', 'DisplayName','estimation data');
    legend;
    ylabel('offset(m)');

    savefig(f, fullfile(temp_folder_path, plots_folder_name,...
                        strcat('SnippetsLDM_', num2str(i), '_', num2str(shiftID), name(1:end-4), '.fig')));
    saveas(f, fullfile(temp_folder_path, plots_folder_name,...
                    strcat('SnippetsLDM_', num2str(i), '_', num2str(shiftID), name(1:end-4), '.png')));
    close(f);

    if (GENERATE_OFFSET_TIME_PLOTS)
        estimationEstimation = input*P';
        f = figure();
        plot(0:dT:dT*(size(estimationEstimation,1)-1),estimationEstimation(:,1)*alfa(1)+mu(1),'LineWidth', 2, 'color', 'k', 'LineStyle','--', 'DisplayName', 'estimated by LDM');    
        hold on;
        xlabel('time(s)');
        grid on;
        plot(0:dT:dT*(size(estimationEstimation,1)-1),output(:,1)*alfa(1)+mu(1), 'LineWidth', 2, 'LineStyle', ':', 'color', 'b', 'DisplayName','validation data');
        legend;
        ylabel('offset(m)');
        xlabel('time(s)');

        savefig(f, fullfile(temp_folder_path, plots_folder_name,...
                            strcat('TimePlotLDM_', num2str(i), name(1:end-4), '.fig')));
        saveas(f, fullfile(temp_folder_path, plots_folder_name,...
                        strcat('TimePlotLDM_', num2str(i), name(1:end-4), '.png')));
        close(f);
    end

end

if (RATIO_OF_TRAIN_VS_TOTAL < 1)
    KPIsum = [mean(KPI(:,2)) max(KPI(:,2)) std(KPI(:,3))];

    f = figure();
    subplot(3,1,1);
    plot(KPI(:,1), 'Marker','x', 'LineWidth', 1.5, 'color', 'b');
    grid on;
    title(strcat('RMS value for', {' '}, num2str(shiftOnOutput), {' '}, 'shift on output'));
    ylabel('offset(m)'); xlabel('epoch');

    subplot(3,1,2);
    plot(KPI(:,2), 'Marker','x', 'LineWidth', 1.5, 'color', 'b');
    grid on;
    title(strcat('NRMS_W value for', {' '}, num2str(shiftOnOutput), {' '}, 'shift on output'));
    xlabel('epoch');

    subplot(3,1,3);
    plot(KPI(:,3), 'Marker','x', 'LineWidth', 1.5, 'color', 'b');
    grid on;
    title(strcat('NRMS_M value for', {' '}, num2str(shiftOnOutput), {' '}, 'shift on output'));
    xlabel('epoch');

    savefig(f, fullfile(temp_folder_path, plots_folder_name,...
                            strcat('KPI_LDM_', name(1:end-4), '.fig')));
    saveas(f, fullfile(temp_folder_path, plots_folder_name,...
                    strcat('KPI_LDM_', name(1:end-4), '.png')));
    close(f);

    save( fullfile(temp_folder_path, plots_folder_name,'KPI_LDM.mat'), 'KPI');
    save( fullfile(temp_folder_path, plots_folder_name,'KPIsum_LDM.mat'), 'KPIsum');
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

function [input, mu, alfa] = normAndCentral(input)
for i=1:size(input,2)
    mu(i) = mean(input(:,i));
    alfa(i) = max(abs(input(:,i)));
    input(:,i) = (input(:,i)-mu(i));
    input(:,i) = input(:,i)/alfa(i);
end
end
