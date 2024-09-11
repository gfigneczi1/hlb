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
OUTPUT_SHIFT = 15:OUTPUT_STEP_DISTANCE:MAXIMUM_PREVIEW_DISTANCE; % preview distance is divided into sections where the output will be predicted
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


%% loop through data sections
% PARAMETERS
p = RATIO_OF_TRAIN_VS_TOTAL; %percentage of evaluation data from entire dataset

for driverID = 1:size(segments.segments,2)
    DRIVER_ID_IF_NOT_MERGED = driverID;
    segment = segments.segments(DRIVER_ID_IF_NOT_MERGED).segment;
    name = segments.segments(DRIVER_ID_IF_NOT_MERGED).name;
    [~, segment_m, indexes] = prepareInputForPlanner(segment);
    dT = mean(diff(segment_m(:, indexes.Relative_time)));
    for shiftID=1:numel(OUTPUT_SHIFT)
        tic;
        dx = segment_m(:, indexes.VelocityX)*dT;
        shiftOnOutput = [1:1:size(segment_m(:,indexes.Relative_time),1)]'+floor(OUTPUT_SHIFT(shiftID)./dx);
        input = [segment_m(:, indexes.OncomingVehicleTimeToPass), ...
                segment_m(:, indexes.VelocityX), ...
                movmean(segment_m(:, indexes.AccelerationX),20), ...
                movmean(segment_m(:, indexes.YawRate),20), ...
                movmean(segment_m(:, indexes.LaneCurvature), 20), ...
                movmean(segment_m(:, indexes.c3), 200)];
    
        variables = ["$t_{pass}$", "$v_x$", "$a_x$", "$\omega$", "$\kappa_{road}$", "$d\kappa$"];
    
        output = zeros(size(input,1),1);
        for shiftIDonOutput=1:size(input,1)
            if (shiftOnOutput(shiftIDonOutput) > size(input,1))
                break;
            else
                output(shiftIDonOutput,1) = -segment_m(shiftOnOutput(shiftIDonOutput), indexes.c0);
            end
        end
        
        outputRaw = output;

    
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
                            strcat('LRM_inputOutputs', name(1:end-4), '_', num2str(i), '.fig')));
        saveas(f, fullfile(temp_folder_path, plots_folder_name,...
                        strcat('LRM_inputOutputs_', name(1:end-4), '_', num2str(i), '.png')));
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
        [output,m_out,s_out] = normAndCentr(output);
        [input,m_in,s_in] = normAndCentr(input);
        
        Din = input;
        Dout = output;
        M = size(Din,1);
        % EVALUATION / VALIDATION DATA SELECTION
        estimationData = 1:1:round(p*M);
        validationData = round(p*M)+1:1:M;
        input_estimation = Din(estimationData,:);
        output_estimation = Dout(estimationData,1);
        input_validation = Din(validationData,:);
        output_validation = Dout(validationData,:);
    
        %% LDM regression
        P = output_estimation'*input_estimation*inv(input_estimation'*input_estimation);
        
        if (RATIO_OF_TRAIN_VS_TOTAL < 1)
            % Evaluation of validation data
            estimation = input_validation*P';
            %% KPI-s
            % RMS calculation
            RMS = sqrt(sum((estimation-output_validation).^2)/size(estimation,1));
            % NRMS - normalization on scale
            W = max(output_validation) - min(output_validation);
            NRMS_W = RMS/W;
            % NRMS - normalization on absolute maximum
            M = max(abs(output_validation));
            NRMS_M = RMS/M;
            KPI(shiftID,1:3) = [RMS NRMS_W NRMS_M];
                    
            fprintf("EVALUTATION of snippet %d:\n", i);
            fprintf("RMS value is: %f\n", RMS);
            fprintf("NRMS value based on range: %f\n", NRMS_W);
            fprintf("NRMS value based on absolute maximum: %f\n", NRMS_M);
        else
            KPI(shiftID,1:3) = [0 0 0];
        end
        
        %% CHECKING TRAIN DATA
        % Evaluation of estimation data
        estimationEstimation = input_estimation*P';
    
        %% KPI-s
        % RMS calculation
        RMSest = sqrt(sum((estimationEstimation-output_estimation).^2)/size(estimationEstimation,1));
        % NRMS - normalization on scale
        W = max(output_estimation) - min(output_estimation);
        NRMS_W_est = RMSest/W;
        % NRMS - normalization on absolute maximum
        M = max(abs(output_estimation));
        NRMS_M_est = RMSest/M;
        
        KPI(shiftID,4:6) = [RMSest NRMS_W_est NRMS_M_est];
        
        %Adding the parameters
        KPI(shiftID,7:6+size(input_estimation,2)) = P;

        fprintf("EVALUTATION of snippet %d:\n", i);
        fprintf("RMS_est value is: %f\n", RMSest);
        fprintf("NRMS_est value based on range: %f\n", NRMS_W_est);
        fprintf("NRMS_est value based on absolute maximum: %f\n", NRMS_M_est);

        f = figure('units','normalized','outerposition',[0 0 1 1]);
        if (RATIO_OF_TRAIN_VS_TOTAL < 1)
            % plot the estimation data
            subplot(2,1,1);
            hold on;
            plot(estimation(:,1),'LineWidth', 2, 'color', 'k', 'LineStyle','--', 'DisplayName', 'estimated by LRM');    
            ylabel('offset');
            grid on;
            plot(output_validation(:,1), 'LineWidth', 2, 'LineStyle', ':', 'color', 'b', 'DisplayName','validation data');
            legend;
            ylabel('offset(m)');
        end

        % plot the estimation data
        subplot(2,1,2);
        hold on;
        plot(estimationEstimation(:,1),'LineWidth', 2, 'color', 'k', 'LineStyle','--', 'DisplayName', 'estimated by LRM');    
        ylabel('offset');
        grid on;
        plot(output_estimation(:,1), 'LineWidth', 2, 'LineStyle', ':', 'color', 'b', 'DisplayName','estimation data');
        legend;
        ylabel('offset(m)');

        savefig(f, fullfile(temp_folder_path, plots_folder_name,...
                            strcat('SnippetsLRM_', num2str(i), '_', num2str(shiftID), name(1:end-4), '.fig')));
        saveas(f, fullfile(temp_folder_path, plots_folder_name,...
                        strcat('SnippetsLRM_', num2str(i), '_', num2str(shiftID), name(1:end-4), '.png')));
        close(f);
        
        if (GENERATE_OFFSET_TIME_PLOTS)
            
            estimationEstimation = input*P';
            estimationEstimation = estimationEstimation*diag(s_out)+m_out;
            
            f = figure();
            f.Position = [100 100 650 550];
            set(f,'defaulttextInterpreter','latex') ;
            set(f, 'defaultAxesTickLabelInterpreter','latex');  
            set(f, 'defaultLegendInterpreter','latex');
            subplot(2,1,1);
            plot(0:dT:dT*(size(estimationEstimation,1)-1),estimationEstimation(:,1),'LineWidth', 2, 'color', 'k', 'LineStyle',':', 'DisplayName', 'Estimated by LRM');    
            hold on;
            xlabel('$t(s)$');
            grid on;
            plot(0:dT:dT*(size(estimationEstimation,1)-1),output(:,1)*s_out+m_out, 'LineWidth', 2, 'LineStyle', ':', 'color', 'b', 'DisplayName','Measured data');
            legend;
            ylabel('$\delta_1(m)$');
            xlabel('$t(s)$');
            
            subplot(2,1,2);
            plot(0:dT:dT*(size(estimationEstimation,1)-1), input(:,5)*s_in(5)+m_in(5), 'LineWidth',2,'color', 'k', 'DisplayName', 'Curvature');
            ylabel('$\delta_1(m)$');
            xlabel('$t(s)$');            
            
            savefig(f, fullfile(temp_folder_path, plots_folder_name,...
                                strcat('TimePlotLRM_', num2str(i), name(1:end-4), '.fig')));
            saveas(f, fullfile(temp_folder_path, plots_folder_name,...
                            strcat('TimePlotLRM_', num2str(i), name(1:end-4), '.png')));
            close(f);
        end
    end
    save( fullfile(temp_folder_path, plots_folder_name,strcat('KPI_LRM_driver_', num2str(driverID), '.mat')), 'KPI');
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

