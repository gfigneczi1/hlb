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
RATIO_OF_TRAIN_VS_TOTAL = 0.5;
GENERATE_OFFSET_TIME_PLOTS = false;
MAXIMUM_PREVIEW_DISTANCE = 150; % from within preview information is extracted, unit: m
OUTPUT_STEP_DISTANCE = 10; % in meters
NUMBER_OF_PREVIEW_INFORMATION = 2; % maximum number
OUTPUT_SHIFT = linspace(15,MAXIMUM_PREVIEW_DISTANCE,10); %10:OUTPUT_STEP_DISTANCE:MAXIMUM_PREVIEW_DISTANCE; % preview distance is divided into sections where the output will be predicted
MERGE_DATA_TABLES = false;
DRIVER_ID_IF_NOT_MERGED = 1;
FILTER_OUTPUT = false;

for i=1:length(segments.segments)-2
    % loop through segment_ms and concatenate each to the previous one
    if (i==1)
        segmentMerged=segments.segments(i).segment;
    else
        segmentMerged=[segmentMerged; segments.segments(i).segment];
    end
end
[~, segmentMerged_m, indexesMerged] = prepareInputForPlanner(segmentMerged);

PARAMS.MAXIMUM_INPUT_SIZE = MAXIMUM_INPUT_SIZE;
PARAMS.FILTER_OUTPUT = FILTER_OUTPUT;
PARAMS.OUTPUT_SHIFT = OUTPUT_SHIFT;
[~, ~, ~, ~, ~, s_out, ~, s_in] = prepareData(segmentMerged_m, indexesMerged, PARAMS);

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

for driverID = 1:size(segments.segments,2)-2
    DRIVER_ID_IF_NOT_MERGED = driverID;
    segment = segments.segments(DRIVER_ID_IF_NOT_MERGED).segment;
    name = segments.segments(DRIVER_ID_IF_NOT_MERGED).name;
    % transforming the struct to array for lower calculation time
    [~, segment_m, indexes] = prepareInputForPlanner(segment);
    dT = mean(diff(segment_m(:, indexes.Relative_time)));

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

        variables = ["$t_{pass}$", "$o_{type}$", "$fo_{type}$", "$v_x$", "$a_x$", "$\omega$", "$\kappa_{road}$", "$d\kappa$"];

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

        % LIMIT DATA IF NEEDED
        % this is done before norm and central
        input = input(1:min(size(input,1),MAXIMUM_INPUT_SIZE),:);
        output = output(1:min(size(input,1),MAXIMUM_INPUT_SIZE), :);

        % NORM AND CENTRAL
        % According to the global cin and sin, cout and sout
        [~, cout, ~] = normalize(output);
        [~, cin, ~] = normalize(input);
        output = (output-cout)./s_out(shiftID);
        input = (input-cin)./s_in;
        

        % REMOVING SPARE POINTS I.E., OUTPUT HIGHER THAN 2, = 2 sigma
        % variance
        input = input(abs(output(:,1))<=2,:);
        output = output(abs(output(:,1))<=2,:);
        MAXIMUM_LENGTH_OF_SNIPPET = length(output);
        
        % SNIPETTING
        for snippetID = 1:min(MAXIMUM_NUMBER_OF_SNIPPET,floor(N/MAXIMUM_LENGTH_OF_SNIPPET))
            snippets{snippetID,1} = input((snippetID-1)*MAXIMUM_LENGTH_OF_SNIPPET+1:snippetID*MAXIMUM_LENGTH_OF_SNIPPET,:);
            snippets{snippetID,2} = output((snippetID-1)*MAXIMUM_LENGTH_OF_SNIPPET+1:snippetID*MAXIMUM_LENGTH_OF_SNIPPET,:);
        end

        % SIMULATION
        for i=1:min(MAXIMUM_NUMBER_OF_SNIPPET,length(snippets))
            Din = snippets{i,1};
            Dout = snippets{i,2};
            M = size(Din,1);
            % EVALUATION / VALIDATION DATA SELECTION
            estimationData = 1:1:round(p*M);
            validationData = round(p*M)+1:1:M;
            input_estimation = Din(estimationData,:);
            output_estimation = Dout(estimationData,1);
            input_validation = Din(validationData,:);
            output_validation = Dout(validationData,:);

            % resampling estimation data
            k = 1;
            input_estimation = input_estimation(1:k:end,:);
            output_estimation = output_estimation(1:k:end,:);
    
            %% Fitting the PP3 model
            PP3_coefficients = polynomialRegression(input_estimation,output_estimation,3);
            estimationEstimation = reFitPolynomial(input_estimation, PP3_coefficients, 3);

            if (RATIO_OF_TRAIN_VS_TOTAL < 1)
                % Evaluation of validation data
                estimation = reFitPolynomial(input_validation, PP3_coefficients, 3);

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
            else
                RMS=0; NRMS_W = 0; NRMS_M = 0; RMS_DEV = 0;
            end
        
            %% KPI-s
            % RMS calculation
            RMSest = sqrt(sum((estimationEstimation-output_estimation).^2)/size(estimationEstimation,1));
            % NRMS - normalization on scale
            W = max(output_estimation) - min(output_estimation);
            NRMS_W_est = RMSest/W;
            % NRMS - normalization on absolute maximum
            M = max(abs(output_estimation));
            NRMS_M_est = RMSest/M;
            
            fprintf("EVALUTATION of snippet %d:\n", i);
            fprintf("RMS_est value is: %f\n", RMSest);
            fprintf("NRMS_est value based on range: %f\n", NRMS_W_est);
            fprintf("NRMS_est value based on absolute maximum: %f\n", NRMS_M_est);
    
            f = figure('units','normalized','outerposition',[0 0 1 1]);
        
            if (RATIO_OF_TRAIN_VS_TOTAL < 1)
                % plot the estimation data
                subplot(2,1,1);
                plot(estimation,'LineWidth', 2, 'color', 'k', 'LineStyle','--', 'DisplayName', 'estimated by GP');    
                ylabel('offset');
                grid on;
                hold on;
                plot(output_validation, 'LineWidth', 2, 'LineStyle', ':', 'color', 'b', 'DisplayName','validation data');
                legend;
                ylabel('offset(m)');
            end
    
            % plot the estimation data
            subplot(2,1,2);
            plot(estimationEstimation,'LineWidth', 2, 'color', 'k', 'LineStyle','--', 'DisplayName', 'estimated by GP');    
            ylabel('offset');
            hold on;
            grid on;
            plot(output_estimation, 'LineWidth', 2, 'LineStyle', ':', 'color', 'b', 'DisplayName','estimation data');
            legend;
            ylabel('offset(m)');
            
            savefig(f, fullfile(temp_folder_path, plots_folder_name,...
                                strcat('Snippets_', num2str(i), '_', num2str(shiftID), '_driver', num2str(driverID), '.fig')));
            saveas(f, fullfile(temp_folder_path, plots_folder_name,...
                            strcat('Snippets_', num2str(i), '_', num2str(shiftID), '_driver', num2str(driverID), '.png')));
            close(f);      
       
            trun = toc;
            KPI{shiftID,i}= [RMS NRMS_W NRMS_M RMSest NRMS_W_est NRMS_M_est PP3_coefficients' trun];
            ETA(shiftID).input_estimation = input_estimation;
            ETA(shiftID).output_estimation = output_estimation;
            ETA(shiftID).normFactors = [cin, s_in, cout,s_out];

            fprintf("time of shift id %d is %f\n", shiftID, trun);        
        end
    
    end % end of shift selection
    save( fullfile(temp_folder_path, plots_folder_name,strcat('KPI_', '_driver_', num2str(driverID), '.mat')), 'KPI');
    save( fullfile(temp_folder_path, plots_folder_name,strcat('ETA_','_driver_', num2str(driverID), '.mat')), 'ETA');
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

function [P_array, delta_0] = trainLDM(input, output, PARAMS)
% LDM requires special input set, therefore a local input_loc is created 
% with necessary data
dT = 0.05;
ds = input(:,2)*dT;
delta_0 = mean(output(:,1));
for i=1:size(input,1)
    S = cumsum(ds(i:end));
    np1 = i+find(S>=PARAMS(1),1);
    np2 = i+find(S>=PARAMS(2),1);
    np3 = i+find(S>=PARAMS(3),1)-1;
    if (isempty(np3))
        break;
    end
    input_loc(i,:) = [mean(input(i:np1,5)) mean(input(np1:np2,5)) mean(input(np2:np3,5))];
    output_loc(i,:) = output([np1, np2, np3],1)-delta_0;
end
P = functional_driverModelLearning(input_loc/0.001, output_loc ,8);
delta_0 = P(19);
P = reshape(P(1:18),3,6);

P_array{1} = P(:,1);
P_array{2} = P(:,2);
P_array{3} = P(:,3);
P_array{4} = P(:,4);
P_array{5} = P(:,5);
P_array{6} = P(:,6);

end

function delta = estimateLDM(input, dT)
    dT = 0.05;
    ds = input(:,2)*dT;
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
    leftCurves = mean(input_loc') >= 2.5e-4;
    rightCurves = mean(input_loc') <= -2.5e-4;
    inputLeft = input_loc; inputLeft(~leftCurves,:) = 0;
    inputRight = input_loc; inputRight(~rightCurves,:) = 0;

    for i=1:length(LDM_params.P_array)/2
        estimationLDM{i} = inputLeft/0.001*LDM_params.P_array{i}+inputRight/0.001*LDM_params.P_array{i+3}+LDM_params.delta_0;
    end
end

function [input, output, inputRaw, outputRaw, c_out, s_out, c_in, s_in] = prepareData(segment_m, indexes, PARAMS)
     % input array
     input = [segment_m(:, indexes.OncomingVehicleTimeToPass), ...
                segment_m(:, indexes.OncomingTrafficType), ...
                segment_m(:, indexes.FrontTrafficType), ...
                segment_m(:, indexes.VelocityX), ...
                movmean(segment_m(:, indexes.AccelerationX),20), ...
                movmean(segment_m(:, indexes.YawRate),20), ...
                movmean(segment_m(:, indexes.LaneCurvature), 20), ...
                movmean(segment_m(:, indexes.c3), 200)];
    
    % output array
    output = zeros(size(input,1),numel(PARAMS.OUTPUT_SHIFT));
    if (PARAMS.FILTER_OUTPUT)
        delta = movmean(-segment_m(:,indexes.c0),180);
    else
        delta = -segment_m(:,indexes.c0);
    end
    for shiftID=1:numel(PARAMS.OUTPUT_SHIFT)
        dT = mean(diff(segment_m(:, indexes.Relative_time)));
        dx = segment_m(:, indexes.VelocityX)*dT;        
        shiftOnOutput = [1:1:size(segment_m(:,indexes.Relative_time),1)]'+floor(PARAMS.OUTPUT_SHIFT(shiftID)./dx);
        for shiftIDonOutput=1:size(input,1)
            if (shiftOnOutput(shiftIDonOutput) > size(input,1))
                break;
            else
                output(shiftIDonOutput,shiftID) = delta(shiftOnOutput(shiftIDonOutput));
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

function coefficients = polynomialRegression(X,y, p)
% formula from here: https://math.stackexchange.com/questions/3155866/multivariate-quadratic-regression
% mathematically: (A'A)^(-1)*A'*y
% where: A[i,:] = [1 xi1 xi2 ... xin xi1^2 xi2^2 ... xin^2 ...xin^p] where
% p is the order of the polynomial
for i = 1:size(X,1)
    % loop through input samples
    x = X(i,:);
    for j=1:p
        A(i,(j-1)*size(X,2)+1:j*size(X,2)) = x.^(j);
    end
end
coefficients = inv(A'*A)*A'*y;
end

function y_ = reFitPolynomial(X, c, p)
y_ = zeros(size(X,1),1);
for i=1:p
    y_ = y_+(c((i-1)*size(X,2)+1:i*size(X,2))'*(X.^i)')';
end
end
