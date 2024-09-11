close all; clear;
load('C:\database\KDP_HLB_GP\Dr008_Dr023_input_GP.mat');
config.root = "./";

MAXIMUM_LENGTH_OF_SNIPPET = 25000;
MAXIMUM_INPUT_SIZE = 25000;
MAXIMUM_NUMBER_OF_SNIPPET = 1;
RATIO_OF_TRAIN_VS_TOTAL = 0.7;
GENERATE_OFFSET_TIME_PLOTS = true;
MAXIMUM_PREVIEW_DISTANCE = 150; % from within preview information is extracted, unit: m
OUTPUT_STEP_DISTANCE = 15; % in meters
MERGE_DATA_TABLES = false;
DRIVER_ID_IF_NOT_MERGED = 5;
SIMPLIFICATION_RATIO = 0;
LDM_NP = [10, 39, 136];

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


%% loop through data sections
% PARAMETERS
p = RATIO_OF_TRAIN_VS_TOTAL; %percentage of evaluation data from entire dataset

for driverID = 1:size(segments.segments,2)
    DRIVER_ID_IF_NOT_MERGED = driverID;
    segment = segments.segments(DRIVER_ID_IF_NOT_MERGED).segment;
    name = segments.segments(DRIVER_ID_IF_NOT_MERGED).name;
    [~, segment_m, indexes] = prepareInputForPlanner(segment);
    dT = mean(diff(segment_m(:, indexes.Relative_time)));

    tic;
    % LDM requires special input set, therefore a local input_loc is created 
    % with necessary data
    dT = 0.05;
    ds = segment_m(:,indexes.VelocityX)*dT;
    outputRaw = movmean(-segment_m(:,indexes.c0),180);
    for i=1:size(segment_m,1)
        S = cumsum(ds(i:end));
        np1 = i+find(S>=LDM_NP(1),1);
        np2 = i+find(S>=LDM_NP(2),1);
        np3 = i+find(S>=LDM_NP(3),1)-1;
        if (isempty(np3))
            break;
        end
        input(i,:) = [mean(segment_m(i:np1,indexes.LaneCurvature)) mean(segment_m(np1:np2,indexes.LaneCurvature)) mean(segment_m(np2:np3,indexes.LaneCurvature))];
        output(i,:) = outputRaw([np1, np2, np3],1);
    end

    variables = ["$\overline{\kappa}_{on}$", "$\overline{\kappa}_{nm}$", "$\overline{\kappa}_{mf}$"];

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
                        strcat('ELDM_inputOutputs', name(1:end-4), '_', num2str(i), '.fig')));
    saveas(f, fullfile(temp_folder_path, plots_folder_name,...
                    strcat('ELDM_inputOutputs_', name(1:end-4), '_', num2str(i), '.png')));
    close(f);

    % SHUFFLE
    N = size(input,1);
    shuffledIndeces = randperm(N);
    input = input(shuffledIndeces,:);
    output = output(shuffledIndeces, :);

    % LIMIT DATA IF NEEDED
    % this is done before norm and central
    input = input(1:min(size(input,1),MAXIMUM_INPUT_SIZE),:);
    output = output(1:min(size(input,1),MAXIMUM_INPUT_SIZE), :);
   
    Din = input;
    Dout = output;
    M = size(Din,1);
    % EVALUATION / VALIDATION DATA SELECTION
    estimationData = 1:1:round(p*M);
    validationData = round(p*M)+1:1:M;
    input_estimation = Din(estimationData,:);
    output_estimation = Dout(estimationData,:);
    input_validation = Din(validationData,:);
    output_validation = Dout(validationData,:);

    %% ELDM regression
    P = functional_driverModelLearning(input_estimation/0.001,output_estimation,8);
    
    if (RATIO_OF_TRAIN_VS_TOTAL < 1)
        % Evaluation of validation data
        leftCurves = mean(input_validation') >= 2.5e-4;
        rightCurves = mean(input_validation') <= -2.5e-4;
        inputLeft = input_validation; inputLeft(~leftCurves,:) = 0;
        inputRight = input_validation; inputRight(~rightCurves,:) = 0;

        estimation = inputLeft/0.001*reshape(P(1:9),3,3)+inputRight/0.001*reshape(P(10:18),3,3)+P(19);
        
        %% KPI-s
        % RMS calculation
        RMS = sqrt(sum((estimation-output_validation).^2)/size(estimation,1));
        % NRMS - normalization on scale
        W = max(output_validation) - min(output_validation);
        NRMS_W = RMS./W;
        % NRMS - normalization on absolute maximum
        M = max(abs(output_validation));
        NRMS_M = RMS./M;
        KPI(1:3,1:3) = [RMS' NRMS_W' NRMS_M'];
                
        fprintf("EVALUTATION of snippet %d:\n", i);
        fprintf("RMS value is: %f\n", RMS);
        fprintf("NRMS value based on range: %f\n", NRMS_W);
        fprintf("NRMS value based on absolute maximum: %f\n", NRMS_M);
    else
        KPI(1:3,1:3) = zeros(3,3);
    end
    
    %% CHECKING TRAIN DATA
    % Evaluation of estimation data
    leftCurves = mean(input_estimation') >= 2.5e-4;
    rightCurves = mean(input_estimation') <= -2.5e-4;
    inputLeft = input_estimation; inputLeft(~leftCurves,:) = 0;
    inputRight = input_estimation; inputRight(~rightCurves,:) = 0;

    estimationEstimation = inputLeft/0.001*reshape(P(1:9),3,3)+inputRight/0.001*reshape(P(10:18),3,3)+P(19);

    %% KPI-s
    % RMS calculation
    RMSest = sqrt(sum((estimationEstimation-output_estimation).^2)/size(estimationEstimation,1));
    % NRMS - normalization on scale
    W = max(output_estimation) - min(output_estimation);
    NRMS_W_est = RMSest./W;
    % NRMS - normalization on absolute maximum
    M = max(abs(output_estimation));
    NRMS_M_est = RMSest./M;
    
    KPI(1:3,4:6) = [RMSest' NRMS_W_est' NRMS_M_est'];
    
    %Adding the parameters
    KPI(1:3,7:13) = [reshape(P(1:18),3,6) ones(3,1)*P(19)];

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
                        strcat('SnippetsELDM_', num2str(i), '_', name(1:end-4), '.fig')));
    saveas(f, fullfile(temp_folder_path, plots_folder_name,...
                    strcat('SnippetsELDM_', num2str(i), '_', name(1:end-4), '.png')));
    close(f);

    save( fullfile(temp_folder_path, plots_folder_name,strcat('KPI_ELDM_driver_', num2str(driverID), '.mat')), 'KPI');
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

function [dataOut,m,s] = normAndCentr(dataIn)
for i=1:size(dataIn,2)
    m(i) = mean(dataIn(:,i));
    s(i) = std(dataIn(:,i));
    dataOut(:,i) = (dataIn(:,i)-m(i))/s(i);
end
end

