function evaluator_fullDriverModel(segments,config)
SIMULATION_MODE = "SYNTHETIC_SIMULATION";
global segment_m indexes globalStartIndex globalStopIndex vehicleState parameters metadata corFine net refFine modelID nrmsLoss nrmsLossManual dts modelMode

% path = [0:1:1000, zeros(1,1001)]';
% targetSpeed = 30; %kph
% vehiclePath = functional_vehicleModelTest(path, targetSpeed);
%     

parameters.P_npDistances = [0.04, 0.116, 0.388];
parameters.numberOfNodePoints = 3;
parameters.drID = 8;
%parameters.P_ELDM = reshape([2.86078399132926;0.884814215672794;-1.90657794718284;-3.09943416608130;-0.665457759838954;2.30236448840005;0.348462602099426;-0.107035325513227;-0.271014703397729;1.07959046302992;-0.775251579323662;-0.252977961446196;-0.822164501814478;1.36747233514778;0.113183483561418;-0.124241139196637;-0.454142531428492;0.293625990988783;-0.000983031283019174;-0.000983031283019174;-0.000983031283019174], 3,7);
parameters.P_ELDM = zeros(3,7);
parameters.P_replanCycle = 10;
parameters.vehicleParameters.wheelBase = 2.7;
parameters.vehicleParameters.r = 0.309725; 
parameters.vehicleParameters.c_alfaf = 8600; 
parameters.vehicleParameters.c_sf = 23000; 
parameters.vehicleParameters.c_alfar = 8600; 
parameters.vehicleParameters.c_sr = 23000;
parameters.vehicleParameters.m = 1519;
parameters.vehicleParameters.Jwheel = 250; 
parameters.vehicleParameters.J = 1818;
parameters.vehicleParameters.A = 1.5; 
parameters.vehicleParameters.c_w = 0; 
parameters.vehicleParameters.rho_air = 1; 
parameters.vehicleParameters.lf = 1; 
parameters.vehicleParameters.lr = 1.5; 

PARAMS.MAXIMUM_INPUT_SIZE = 15000; % before snipetting and normalization
PARAMS.KERNEL_TYPE ="{'covPPard',3}"; % "{'covSum', {'covLINard', {'covPPard',3}}}"; %"{'covSum', {'covSEard', {'covPPard',3}}}";
PARAMS.RATIO_OF_TRAIN_VS_TOTAL = 0.7;
PARAMS.MAXIMUM_PREVIEW_DISTANCE = 150; % from within preview information is extracted, unit: m
PARAMS.OUTPUT_STEP_DISTANCE = 10; % in meters
PARAMS.NUMBER_OF_PREVIEW_INFORMATION = 2; % maximum number
PARAMS.OUTPUT_SHIFT = linspace(15,PARAMS.MAXIMUM_PREVIEW_DISTANCE,10); %10:OUTPUT_STEP_DISTANCE:MAXIMUM_PREVIEW_DISTANCE; % preview distance is divided into sections where the output will be predicted
PARAMS.EPOCH_CALCULATION = true;
PARAMS.EPOCH_NUMBER = 10;
PARAMS.GREEDY_REDUCTION = true;
PARAMS.LDM_NP = [10, 39, 136];
PARAMS.DriverID = 4;
PARAMS.FILTER_OUTPUT = false;

parameters.PARAMS = PARAMS;

switch SIMULATION_MODE
    case "SYNTHETIC_SIMULATION"
        colors = ["r", "b", "c", "g", "k", "m", "r--", "b--"];
        for segmentID = 1:length(segments.segments)
            segment = segments.segments(segmentID).segment;
            name = segments.segments(segmentID).name;
            TC = name(end-5:end-4);
            [~, segment_m, indexes] = prepareInputForPlanner(segment);
            modelMode = "kinematic";
            modelID = "gp";
            
            for driverID=1:8 % currently static as driver parameter read happens inside subfunction!
                parameters.PARAMS.DriverID = driverID;
                %% GPM          
                [estimation, deviation, inputRaw] = gpGenerateEstimate(segment_m, indexes);
                segment_m(:,end+1:end+10) = estimation;
                indexes.GP_1 = size(segment_m,2)-9;
                indexes.GP_10 = size(segment_m,2);
                segment_m(:,end+1:end+10) = deviation;
                indexes.errGP_1 = size(segment_m,2)-9;
                indexes.errGP_10 = size(segment_m,2);
                
                globalStartIndex = 2; % minimum is 2, otherwise it fails!
                globalStopIndex = size(segment_m,1);
        
                [path, ~, ~, ~, intentionPath] = pathGeneration();
                corFine = corridorGeneration(0.5);
                refFine = referenceGeneration(0.5);
        
                offsetPath = offsetCalculation(corFine, path);
                offsetRef = offsetCalculation(corFine, refFine);
        
                offsets{driverID}.offset = offsetPath; offsets{driverID}.name = strcat('GPM-Dr', num2str(driverID)); offsets{driverID}.X = path(:,1); offsets{driverID}.marker = colors(driverID);
            end
            f = plot_offsetPlotsExtended(segment_m,indexes,offsets,[], inputRaw);
            temp_folder_path = "../../_temp";
            plots_folder_name = "plots";
            savefig(f, fullfile(temp_folder_path, plots_folder_name,...
                                strcat('Resimulate_plots_syntheticData_TC', TC, '.fig')));
            saveas(f, fullfile(temp_folder_path, plots_folder_name,...
                            strcat('Resimulate_plots_syntheticData_TC', TC, '.png')));
            close(f);
        end

    case "SINGLE_SIMULATION"
        load('C:\database\KDP_HLB_GP\Dr008_Dr023_input_GP.mat');
        segment = segments.segments(PARAMS.DriverID).segment;
        name = segments.segments(PARAMS.DriverID).name;
        %segment = segments.segments(1).segment;
        [~, segment_m, indexes] = prepareInputForPlanner(segment);
        segment_m(abs(segment_m(:,indexes.LaneOrientation))>0.05,indexes.LaneOrientation) = 0;

        [~, ~, inputRaw, outputRaw, parameters.PARAMS.c_out, parameters.PARAMS.s_out, parameters.PARAMS.c, parameters.PARAMS.s] = prepareData(segment_m, indexes, parameters.PARAMS);
        
        modelMode = "kinematic";

        %% PP3 model
        modelID = "pp3";
        estimationPP3 = pp3GenerateEstimate(segment_m, indexes);
        segment_m(:,end+1:end+10) = estimationPP3;
        indexes.PP3_1 = size(segment_m,2)-9;
        indexes.PP3_10 = size(segment_m,2);

        globalStartIndex = 2; % minimum is 2, otherwise it fails!
        globalStopIndex = size(segment_m,1);

        [path, U, dY, ~, intentionPath] = pathGeneration();
        corFine = corridorGeneration(0.5);
        refFine = referenceGeneration(0.5);

        offsetPath = offsetCalculation(corFine, path);
        offsetRef = offsetCalculation(corFine, refFine);

        offsets{6}.offset = offsetPath; offsets{6}.name = 'PP3'; offsets{6}.X = path(:,1); offsets{6}.marker = 'm';


        %% GPM
        modelID = "gp";
                       
        if (modelID=="gp")
            [estimation, deviation] = gpGenerateEstimate(segment_m, indexes);
            segment_m(:,end+1:end+10) = estimation;
            indexes.GP_1 = size(segment_m,2)-9;
            indexes.GP_10 = size(segment_m,2);
            segment_m(:,end+1:end+10) = deviation;
            indexes.errGP_1 = size(segment_m,2)-9;
            indexes.errGP_10 = size(segment_m,2);
        end
        
        globalStartIndex = 2; % minimum is 2, otherwise it fails!
        globalStopIndex = size(segment_m,1);

        [path, U, dY, ~, intentionPath] = pathGeneration();
        corFine = corridorGeneration(0.5);
        refFine = referenceGeneration(0.5);

        offsetPath = offsetCalculation(corFine, path);
        offsetRef = offsetCalculation(corFine, refFine);

        offsets{1}.offset = offsetPath; offsets{1}.name = 'GPM'; offsets{1}.X = path(:,1); offsets{1}.marker = 'r';

        %% LRM
        modelID = "lrm";
        if (modelID=="lrm")
            outputOffsetted = outputRaw-mean(mean(outputRaw));
            for shiftID=1:numel(PARAMS.OUTPUT_SHIFT)                
                P_LRM = outputOffsetted(:,shiftID)'*inputRaw*inv(inputRaw'*inputRaw);
                segment_m(:,end+1) = inputRaw*P_LRM'+mean(mean(outputRaw));
                if (shiftID==1)
                    indexes.LRM_1 = size(segment_m,2);
                elseif (shiftID==numel(PARAMS.OUTPUT_SHIFT))
                    indexes.LRM_10 = size(segment_m,2);
                end
            end
        end
        [path, U, dY, ~, intentionPath] = pathGeneration();
        offsetPath = offsetCalculation(corFine, path);
        offsets{2}.offset = offsetPath; offsets{2}.name = 'LRM'; offsets{2}.X = path(:,1); offsets{2}.marker = 'b';
      
        %% Reference
        modelID = "groundTruth";
        path = pathGeneration();
        offsetPath = offsetCalculation(corFine, path);
        offsets{4}.offset = offsetPath; offsets{4}.name = 'Reference'; offsets{4}.X = path(:,1); offsets{4}.marker = 'k--';

        %% Filtered reference
        modelID = "groundTruthFiltered";
        path = pathGeneration();
        offsetPath = offsetCalculation(corFine, path);
        offsets{5}.offset = offsetPath; offsets{5}.name = 'Reference filtered'; offsets{5}.X = path(:,1); offsets{5}.marker = 'k:';

        %% ELDM
        % read-in parameters 
        data = load("C:\database\LDM_PARAMS\parametersAllDrivers.mat");
        for j=1:length(data.parameters_out)
            if (data.parameters_out(j).drID == parameters.drID)
                parameters.P_ELDM = functional_driverModelLearning(data.parameters_out(j).U', data.parameters_out(j).dY' ,8);
                break;
            end
        end
        parameters.P_ELDM = [reshape(parameters.P_ELDM(1:18),3,6) ones(3,1)*parameters.P_ELDM(19)];
        modelID = "eldm";
        path = pathGeneration();
        offsetPath = offsetCalculation(corFine, path);
        offsets{3}.offset = offsetPath; offsets{3}.name = 'ELDM'; offsets{3}.X = path(:,1); offsets{3}.marker = 'g';

        plot_offsetPlots(segment_m,indexes,offsets, PARAMS.DriverID);

    case "DATA_GENERATION"
        for i=1:length(segments.segments)
            segment = segments.segments(i).segment;
            [~, segment_m, indexes] = prepareInputForPlanner(segment);
            segment_m(:,indexes.c0) = movmean(segment_m(:,indexes.c0),180);
            segment_m(:,indexes.LaneCurvature) = movmean(segment_m(:,indexes.LaneCurvature),180);

            npDistances(1) = parameters.P_npDistances(1)*250;
            for np=2:length(parameters.P_npDistances)
                npDistances(np) = npDistances(np-1)+parameters.P_npDistances(np)*250;
            end
            U = []; dY=[];
            for j=1:parameters.P_replanCycle:size(segment_m,1)
                % loop through data and calculate local distance at all
                % points
                egoPos = [segment_m(j,indexes.X_abs), segment_m(j,indexes.Y_abs)];
                locDistances = (sum((([segment_m(j:end,indexes.X_abs) segment_m(j:end,indexes.Y_abs)]-egoPos).^2)')).^0.5;
                locIndeces = []; realNpDistances=[];
                for np=1:length(npDistances)
                    if (isempty(find(locDistances>npDistances(np),1)))
                        break;
                    end
                    locIndeces(np) = find(locDistances>npDistances(np),1);
                    realNpDistances(np) = locDistances(locIndeces(np));
                end
                if (length(realNpDistances)~=length(npDistances))
                    break;
                elseif (all(abs(realNpDistances-npDistances)<3))
                    % node points are kind-of-good
                    for np=1:length(npDistances)
                        if (np==1)
                            u(np,1) = mean(segment_m(j:j+locIndeces(np)-1,indexes.LaneCurvature));
                        else
                            u(np,1) = mean(segment_m(j+locIndeces(np-1):j+locIndeces(np)-1,indexes.LaneCurvature));
                        end
                        dy(np,1) = -segment_m(j+locIndeces(np)-1,indexes.c0);
                    end
                    U = [U u]; dY = [dY dy];
                end
            end
            plot(U(1,:), dY(1,:), 'bo');
            parameters_out(i).U = U/0.001;
            parameters_out(i).dY = dY;
            parameters_out(i).drID = str2num(segments.segments(i).name(3:5));
         end
         save(fullfile(config.root,"plots", 'parametersAllDrivers.mat'), 'parameters_out');
    case "PARAMETER_LEARNING"
        for i=1:length(segments.segments)
            segment = segments.segments(i).segment;
            [~, segment_m, indexes] = prepareInputForPlanner(segment);
             
            globalStartIndex = 2; % minimum is 2, otherwise it fails!
            globalStopIndex = size(segment_m,1);
             
            % Definition of metadata
            metadata.pathValidity = 1;
            
            corFine = corridorGeneration(0.5);
            refFine = referenceGeneration(0.5);
            
            modelID = "groundTruth";
            modelMode = "kinematic";
        
            [path, U, dY, ~, intentionPath] = pathGeneration();
            parameters.P_ELDM = reshape(functional_driverModelLearning(U, dY ,8), 3,7);
        
            parameters_out(i).U = U/0.001;
            parameters_out(i).dY = dY;
            parameters_out(i).drID = str2num(segments.segments(i).name(3:5));
         end
         save(fullfile(config.root,"plots", 'parametersAllDrivers.mat'), 'parameters_out');
    case "PARAMETER_VALIDATION_FROM_FILE"
        modelID = "eldm";
        modelMode = "kinematic";
        segment = segments.segments(1).segment;
        [~, segment_m, indexes] = prepareInputForPlanner(segment);
        
        globalStartIndex =2; % minimum is 2, otherwise it fails!
        globalStopIndex = size(segment_m,1);
        
        corFine = corridorGeneration(0.5);
        refFine = referenceGeneration(0.5);
        % read-in parameters 
        data = load("C:\database\LDM_PARAMS\parametersAllDrivers.mat");
        k=1;
        for j=1:length(data.parameters_out)
            if (data.parameters_out(j).drID ~= 8 || data.parameters_out(j).drID ~= 9)
                parameterRawData(:,k) = functional_driverModelLearning(data.parameters_out(j).U', data.parameters_out(j).dY' ,7);
                dataOut.data(k).U = data.parameters_out(j).U;
                dataOut.data(k).dY = data.parameters_out(j).dY;
                dataOut.data(k).P_GT = parameterRawData(:,k);
                dataOut.data(k).drID = data.parameters_out(j).drID;
                k = k+1;
            end
        end
        % clustering
        K = 3;
        rng(1); %doc: https://www.mathworks.com/help/matlab/math/controlling-random-number-generation.html
        normalizedParams = parameterRawData(1:18,:)/max(max(abs(parameterRawData(1:18,:))));
        normalizedParams = [normalizedParams; parameterRawData(19,:)/max(abs(parameterRawData(19,:)))];

        [driverClusters,~,sumD] = kmeans(normalizedParams',K, 'Distance','sqeuclidean', 'MaxIter', 1000, 'Replicates',10); 
        dataOut.driverClusters = driverClusters;
        f = figure(5);
        set(f,'defaulttextInterpreter','latex') ;
        set(f, 'defaultAxesTickLabelInterpreter','latex');  
        set(f, 'defaultLegendInterpreter','latex');
        set(f, "Position", [100,100,505,310]);
        silhouette(normalizedParams',driverClusters);
        s = silhouette(normalizedParams',driverClusters);
        xline(mean(s),'LineWidth',2); grid on;
        title(strcat("Clustering performance K = ", {' '}, num2str(K),", Manhatten norm"));
        legend("Elements", "Mean");
        set(gca,'FontSize', 14);
        
        %% batch parameter simulation
        clusterIDs = unique(driverClusters);
        f = figure();
        set(f,'defaulttextInterpreter','latex') ;
        set(f, 'defaultAxesTickLabelInterpreter','latex');  
        set(f, 'defaultLegendInterpreter','latex');
        f.Position = [100 100 650 850];

        for k=1:length(clusterIDs)
            subplot(length(clusterIDs),1,k);
            noDriversInCluster = 0;
            P_averaged = [];
            numberOfDrivers = size(parameterRawData,2);            
            for j=1:numberOfDrivers
                drID = data.parameters_out(j).drID;
                driverParam = parameterRawData(:,j);
                if (driverClusters(j) == clusterIDs(k))
                   noDriversInCluster = noDriversInCluster + 1;
                   P = reshape(driverParam(1:18),3,6);
                   parameters.P_ELDM = [P [driverParam(end); driverParam(end); driverParam(end)]];
                   if (isempty(P_averaged))
                       P_averaged = parameters.P_ELDM;
                       Ptemp = reshape(parameters.P_ELDM,21,1);
                       P_all = Ptemp(1:19);
                   else
                       P_averaged = 1/noDriversInCluster * ((noDriversInCluster-1)*P_averaged + parameters.P_ELDM);
                       Ptemp = reshape(parameters.P_ELDM,21,1);
                       P_all = [P_all Ptemp(1:19)];
                   end
                   [path, Udrivers{j}, dYdrivers{j}, plannedPath, ~, vehicleStateMemory] = pathGeneration();
                   offsetPath = offsetCalculation(corFine, path);
                   plot(path(:,1), offsetPath, 'DisplayName', strcat("Dr",num2str(drID))); hold on; grid on;
                   title(strcat('Cluster',{' '}, num2str(clusterIDs(k))));
                end
            end
            P_median = median(P_all');
            Ptemp = reshape([P_median P_median(end) P_median(end)], 3,7);
            parameters.P_ELDM = P_averaged; %Ptemp;
            Ptemp = [P_median P_median(end) P_median(end)]; %;
            Pcentroids{k} = reshape(P_averaged,21,1); %Ptemp(1:19);
            disp(P_averaged);
           [path, Uarray{k}, dYarray{k}, plannedPath] = pathGeneration();
           offsetPath = offsetCalculation(corFine, path);
           plot(path(:,1), offsetPath, 'DisplayName', 'Averaged', 'LineWidth', 2, 'color', 'k');
           xlabel('UTF-X(m)'); ylabel('offset(m)');
           dataSaveForVisualization(config, Uarray{k}, dYarray{k}, path, num2str(clusterIDs(k)), segment_m, indexes,offsetPath, plannedPath, P_averaged);
           yline(0,"HandleVisibility","off");
    
            legend("Orientation","horizontal","Location","southoutside", "NumColumns",5);
            ylim([-1.25, 1.25]);
            xlim([min(path(~isnan(path(:,1)),1)), max(path(~isnan(path(:,1)),1))]);
            yticks([-1,-0.5,0,0.5,1]);
            set(gca,'TickLabelInterpreter','latex');
            set(gca,'FontSize', 12);
        end
end

dataOut.Pcentroids = Pcentroids;
save(fullfile(config.root,"plots", 'dataOut.mat'), 'dataOut');

f = figure(3);
set(f,'defaulttextInterpreter','latex') ;
set(f, 'defaultAxesTickLabelInterpreter','latex');  
set(f, 'defaultLegendInterpreter','latex');
f.Position = [100 100 505 310];
ulims = linspace(-5,5, 100);
markers = ['o', 'x', '*'];
colors = ['b', 'r', 'g'];
for i=1:length(Pcentroids)
    plot(Uarray{i}(1,:),dYarray{i}(1,:), 'color', colors(i), 'Marker', markers(i), 'LineStyle', 'none', 'HandleVisibility','off'); hold on;
    plot(Uarray{i}(4,:),dYarray{i}(1,:), 'color', colors(i), 'Marker', markers(i), 'LineStyle', 'none', 'HandleVisibility','off');
    % fit curve
    cLeft = polyfit(Uarray{i}(1,Uarray{i}(1,:)>0.25),dYarray{i}(1,Uarray{i}(1,:)>0.25),1);
    cRight = polyfit(Uarray{i}(4,Uarray{i}(4,:)<0.25),dYarray{i}(1,Uarray{i}(4,:)<0.25),1);
    plot(ulims(51:100),cLeft(2)+cLeft(1)*ulims(51:100), 'LineWidth',2, 'color', colors(i), 'DisplayName',strcat("Cluster",{' '},num2str(i))); 
    plot(ulims(1:50),cRight(2)+cRight(1)*ulims(1:50), 'LineWidth',2, 'color', colors(i), 'HandleVisibility','off');
end
grid on;
legend('Orientation','horizontal','Location','north');
title("Correlation plot of clusters");
ylim([-1.25,1.25]);
xlim([-5,5]);
xlabel("Normalized curvature to 0.001 (1/m)");
ylabel("Lane offset (m)");
set(gca,'FontSize',14);


    %% loop simulation option for different number of node points
    for i=1:10
        parameters.P_npDistances = ones(1,i)*(140/i)/250;
    
        %parameters.P_npDistances = 135/250; %ones(1,10)*15/250; % [0, 29, 96]/250;
        parameters.P_LDM = reshape([2.75127773046959;-0.0398873783094187;-1.17733327762608;-2.94882588015453;0.460440066139664;1.33967800733396;0.327136469963189;-0.296796325338061;-0.0234709186554829;0;0;0;0;0;0;0;0;0;0;0;0], 3,7);
        parameters.P_LDM = parameters.P_LDM(:,1:3);
        parameters.P_LDM = zeros(length(parameters.P_npDistances));
        parameters.P_ELDM = reshape([2.86078399132926;0.884814215672794;-1.90657794718284;-3.09943416608130;-0.665457759838954;2.30236448840005;0.348462602099426;-0.107035325513227;-0.271014703397729;1.07959046302992;-0.775251579323662;-0.252977961446196;-0.822164501814478;1.36747233514778;0.113183483561418;-0.124241139196637;-0.454142531428492;0.293625990988783;-0.000983031283019174;-0.000983031283019174;-0.000983031283019174], 3,7);
        
        %% Definition of metadata
        metadata.pathValidity = 1;
        
        corFine = corridorGeneration(0.5);
        refFine = referenceGeneration(0.5);
         % Init vehicle state
        vehicleState = initVehicleState(segment_m, indexes, globalStartIndex);
        path = pathGeneration();
        [orientation, curvature] = calcPathGeometry(path);
           
        offsetPath = offsetCalculation(corFine, path);
        
        curves = abs(movmean(curvature,200)) > 2.5e-4;
        
        % loss calculation of resulting path
        indeces = (curves)&(~isnan(offsetPath'));
        loss = median(abs(offsetPath(indeces==1)))
        mean(dts)
        losses(i) = loss;
        times(i) = mean(dts);

    end
    
    disp(losses);
    disp(times);    
end

function dataSaveForVisualization(config, GT_U, GT_Y, path, token, segment_m, indexes, offsetPath, plannedPath, Poptim)
    replan_array = ones(size(path,1),1);
    segment.X_abs = path(:,1);
    segment.Y_abs = path(:,2);
    segment.q_T0 = segment_m(:,indexes.Relative_time);
    segment.c2 = segment_m(:,indexes.LaneCurvature);
    x = linspace(1,size(GT_U,2),size(GT_U,2));
    xFine = linspace(1,size(GT_U,2),size(path,1));
    trajFull = plannedPath;
    corFull = [segment_m(:,indexes.corrX) segment_m(:, indexes.corrY)];
    refFull = path;
    
    for i=1:size(GT_U,1)
        GT_U_array(i,:) = zoh(x,GT_U(i,:),xFine);
    end
    for i=1:size(GT_Y,1)
        GT_Y_array(i,:) = zoh(x,GT_Y(i,:), xFine);
    end
    
    save(fullfile(config.root, strcat(token,'forVisualization.mat')));
end

function offset = offsetCalculation(corFine, path)
    % fine gridding the path
    offset(1:2) = nan;
    for i=3:size(path,1)
        % Transformation matrix
        theta = atan2((path(i,2)-path(i-1,2)),(path(i,1)-path(i-1,1)));
        % searching the closest corridor point on fine grid
        distances = ((corFine(:,1)-path(i,1)).^2+(corFine(:,2)-path(i,2)).^2).^0.5;
        closestIndex = find(distances == min(distances),1);
        clear distances
        if (~isempty(closestIndex) && closestIndex(1,1)<=size(corFine,1)-2)
            % Transformation matrix
            T = [cos(theta) sin(theta); -sin(theta) cos(theta)];
            corFineLocal = (corFine(closestIndex(1,1)-2:closestIndex(1,1)+2,:) - path(i,:))*T';
            x_fine = min(corFineLocal(:,1)):0.005:max(corFineLocal(:,1));
            corFineFineLocal = spline(corFineLocal(:,1), corFineLocal(:,2), x_fine);
            % search for the minimum X distance point
            offset(i) = -corFineFineLocal(find(abs(x_fine)==min(abs(x_fine)),1));
        else
            offset(i) = nan;
        end
    end
end

function corFine = corridorGeneration(step)
    global segment_m indexes
    if (numel(find((diff(segment_m(:,indexes.corrX)))<0))>0)
        if(numel(find((diff(segment_m(:,indexes.corrX)))>0))==0)
            % monotonous decreasing
            x_fine = segment_m(end,indexes.corrX):step:segment_m(1,indexes.corrX);
            corrX = segment_m(:,indexes.corrX);
            corrY = segment_m(:,indexes.corrY);
        elseif(numel(find((diff(segment_m(:,indexes.corrX)))>0))<5)
            x_fine = segment_m(end,indexes.corrX):step:segment_m(1,indexes.corrX);
            corrX = segment_m(find((diff(segment_m(:,indexes.corrX)))<0),indexes.corrX);
            corrY = segment_m(find((diff(segment_m(:,indexes.corrX)))<0),indexes.corrY);
        else
            x_fine=[];
            corrX = segment_m(:,indexes.corrX);
            corrY = segment_m(:,indexes.corrY);
        end
    elseif (numel(find((diff(segment_m(:,indexes.corrX)))>0))>0)
        if(numel(find((diff(segment_m(:,indexes.corrX)))<0))==0)
            % monotonous increasing
            x_fine = segment_m(1,indexes.corrX):step:segment_m(end,indexes.corrX);
            corrX = segment_m(:,indexes.corrX);
            corrY = segment_m(:,indexes.corrY);
        elseif (numel(find((diff(segment_m(:,indexes.corrX)))<0))<5)
            x_fine = segment_m(1,indexes.corrX):step:segment_m(end,indexes.corrX);
            corrX = segment_m(find((diff(segment_m(:,indexes.corrX)))>0),indexes.corrX);
            corrY = segment_m(find((diff(segment_m(:,indexes.corrX)))>0),indexes.corrY);
        else
            x_fine = [];
            corrX = segment_m(:,indexes.corrX);
            corrY = segment_m(:,indexes.corrY);
        end
    else
        %non-monotonous
        x_fine = [];
    end
    
    if (~isempty(x_fine))
        y_fine = spline(corrX, corrY, x_fine');
        corFine = [x_fine' y_fine];
    else
        corFine = [];
    end
end

function refFine = referenceGeneration(step)
    global segment_m indexes
    if (numel(find((diff(segment_m(:,indexes.X_abs)))>0))==0 && ...
        numel(find((diff(segment_m(:,indexes.X_abs)))<0))>0)
        % monotonous decreasing
        x_fine = segment_m(end,indexes.X_abs):step:segment_m(1,indexes.X_abs);
    elseif (numel(find((diff(segment_m(:,indexes.X_abs)))>0))>0 && ...
        numel(find((diff(segment_m(:,indexes.X_abs)))<0))==0)
        % monotonous increasing
        x_fine = segment_m(1,indexes.X_abs):step:segment_m(end,indexes.X_abs);
    else
        %non-monotonous
        x_fine = [];
    end
    
    if (~isempty(x_fine))
        y_fine = spline(segment_m(:,indexes.X_abs), segment_m(:,indexes.Y_abs), x_fine');
        refFine = [x_fine' y_fine];
    else
        refFine = [];
    end
end

function [path, U, dY, plannedPath, intentionPath, vehicleStateMemory] = pathGeneration()
    global segment_m indexes globalStartIndex globalStopIndex vehicleState parameters metadata net dt dts modelMode x_k1
    % Entire model is calculated in the global frame!
    vehicleState = initVehicleState(segment_m, indexes, globalStartIndex, modelMode);
    replan = true;
    scenario = []; % initializing the scenario
    path(1,1) = vehicleState.X;
    path(1,2) = vehicleState.Y;
    replanCounter = 0;
    U = []; dY = []; plannedPath = []; dts = []; intentionPath = [];
    priorPath = zeros(1,2);
    u = zeros(3,1); dy = zeros(3,1); 
    x_k1 = zeros(4,1);
    for i=globalStartIndex:globalStopIndex
        
        % calculate time step
        dT = segment_m(i, indexes.Relative_time) - segment_m(i-1, indexes.Relative_time);
        
        % if there is a valid scenario ongoing, let's check, if there
        % is enough remaining part of it!
        if(~isempty(scenario))
            scenarioFinished = scenarioFinishChecker(posteriorPath, vehicleState);
            if (scenarioFinished)
                replan = 1;
            end
        end
        
        if (replan)
            replanCounter = 0;
            % get the nearest index
            nearestIndex = getNearestIndex(segment_m(:,indexes.X_abs),segment_m(:,indexes.Y_abs), [vehicleState.X vehicleState.Y]);        
            % cut scenario and corresponding information
            % -- parameter: info from past and future = window
            window = [0 160];
            scenario = cutScenario(segment_m, indexes, nearestIndex, vehicleState, window);
            % Planner: incl. planning driver model, intended behavior
            % -- input: the scenario input data from the cutting
            % -- ouput: trajectory points
            if(~isempty(scenario))
                % valid cutting happened, planning can happen too
                try
                    [priorPath, u, dy] = planner (scenario, indexes, path, parameters, net);
                    replan = 0;
                    U = [U u]; dY = [dY dy];
                catch
                    % input data is corrupt and planner failed
                    scenario = [];
                    replan = 1;
                end
                
                dts = [dts; dt];
            end
        end
        % from this point on, only calculate if scenario is valid!!
        if (~isempty(scenario) && dT < 0.1)
            % curve policy: a model which takes the input data and outputs the
            % modified points
            % - cuts the same subsegment as done for the curve
            % get the nearest index again for the moved vehicle
            %nearestIndex = getNearestIndex(segment_m, indexes, vehicleState);
            % do the cut...
            scenarioCurvePolicy = scenario; %cutScenario(segment_m, indexes, nearestIndex, vehicleState, window);
            % calculate the posterior path based on the curve policy
            posteriorPath = priorPath; %curvePolicy (scenario, indexes, priorPath, vehicleState);
            
            % controller
            % - outputs the desired steering angle
            if (modelMode == "dynamic")
                vehicleState = loadModel(vehicleState);
                vehicleState = speedController(vehicleState);
            end

            try
            vehicleState.steeringAngle = controller(posteriorPath, vehicleState); 
            catch
                i;
            end

            % vehicle model
            vehicleState = vehicleModel(vehicleState, dT, modelMode, parameters.vehicleParameters);
            
            % metdata update
            metadata.pathValidity(i) = 1;
            replanCounter = replanCounter + 1;
            if (replanCounter == parameters.P_replanCycle)
                replan = 1;
                replanCounter = 0;
            end
        else
            scenario = [];
            % When there is a corruption in the original data, then set the
            % vehicle position to the closest position in the measurement
            [vehicleState, nearestIndex] = setVehiclePositionToNearest(segment_m, indexes, vehicleState);
            if (isempty(nearestIndex))
                disp('Measurement is corrupt, stopping simulation');
                break;
            end
            replan = 1;
            metadata.pathValidity(i) = 0;
            replanCounter = 0;
            priorPath = path(end,:);
            posteriorPath = priorPath;
        end
        path = [path; [vehicleState.X vehicleState.Y]];
        intentionPath = [intentionPath; priorPath(1,:)];
        plannedPath = [plannedPath; posteriorPath(1,:)];
        vehicleStateMemory{i} = vehicleState;
    end
end

function [path, U, dY, plannedPath, intentionPath, vehicleStateMemory] = pathGenerationLite()
    global segment_m indexes globalStartIndex globalStopIndex vehicleState parameters metadata net dt dts modelMode x_k1
    % Entire model is calculated in the global frame!
    vehicleState = initVehicleState(segment_m, indexes, globalStartIndex, modelMode);
    vehicleStateCheck = initVehicleState(segment_m, indexes, globalStartIndex, "dynamicSimplified");
    replan = 10;
    scenario = []; % initializing the scenario
    path(1,1) = vehicleState.X;
    path(1,2) = vehicleState.Y;
    replanCounter = 0;
    U = []; dY = []; plannedPath = []; dts = []; intentionPath = [];
    u = zeros(3,1); dy = zeros(3,1); 
    coefficients_k1 = [];
    coefficients = coefficients_k1;
    displacementVector = zeros(2,1);
    x_k1 = zeros(5,1); x_k1(3) = segment_m(1,indexes.VelocityX);

    for i=globalStartIndex:globalStopIndex
        % calculate time step
        dT = segment_m(i, indexes.Relative_time) - segment_m(i-1, indexes.Relative_time);
        
        % if there is a valid scenario ongoing, let's check, if there
        % is enough remaining part of it!
        if(~isempty(scenario))
            scenarioFinished = coefficients(end).sectionBorders(2) < 10; % if we are within a certain distance to the last node point, re-init planner
            if (scenarioFinished)
                replan = 1;
            end
        end
        
        if (replan)
            replanCounter = 0;
            % get the nearest index
            nearestIndex = getNearestIndex(segment_m(:,indexes.X_abs),segment_m(:,indexes.Y_abs), [vehicleState.X vehicleState.Y]);        
            % cut scenario and corresponding information
            % -- parameter: info from past and future = window
            window = [0 150];
            scenario = cutScenario(segment_m, indexes, nearestIndex, vehicleState, window);
            % Planner: incl. planning driver model, intended behavior
            % -- input: the scenario input data from the cutting
            % -- ouput: trajectory points
            if(~isempty(scenario))
                % valid cutting happened, planning can happen too
                % first, transform the corridor coefficients to the
                % simulated ego frame
                deltaPose = [scenario(1,indexes.X_abs)-vehicleState.X, ...
                    scenario(1, indexes.Y_abs)-vehicleState.Y, ...
                    scenario(1,indexes.theta_calc)-vehicleState.theta];
                [c0, c1] = transformVideoData(scenario(1,indexes.c0),scenario(1,indexes.LaneOrientation), ...
                    [scenario(1,indexes.X_abs), scenario(1,indexes.Y_abs), scenario(1,indexes.theta_calc)], ...
                    [vehicleState.X, vehicleState.Y, vehicleState.theta]);
                
                [coefficients, u, dy, X,Y,theta] = eldmLite (scenario, indexes, vehicleState, coefficients, parameters.P_npDistances, parameters.P_ELDM, c0, c1);
                replan = 0;
                U = [U u]; dY = [dY dy];
                origo = [vehicleState.X, vehicleState.Y, vehicleState.theta]'; % in [m m rad]
                dts = [dts; dt];
            end
        else
            % transforming the coefficients to the ego frame
            [coefficients, ~] = transformCoefficients(X,Y,theta, origo, [vehicleState.X vehicleState.Y vehicleState.theta]);
        end
        % from this point on, only calculate if scenario is valid!!
        if (~isempty(scenario) && dT < 0.1)
            % controller
            % - outputs the desired steering angle
            if (modelMode == "dynamic")
                vehicleState = loadModel(vehicleState);
                vehicleState = speedController(vehicleState);
            end
            vehicleState.steeringAngle = controllerLite(coefficients, vehicleState, dT); 
            vehicleStateCheck.steeringAngle = vehicleState.steeringAngle;
            % vehicle model
            vehicleState = vehicleModel(vehicleState, dT, modelMode, parameters.vehicleParameters);            
%             vehicleStateCheck = vehicleModel(vehicleStateCheck, dT, "dynamicSimplified");
%             vehicleStateCheck = checkModelEquality(vehicleState,vehicleStateCheck);
%             
%             subplot(2,1,1);
%             plot(vehicleStateCheck.X,vehicleStateCheck.Y, 'bo', ...
%                 vehicleState.X, vehicleState.Y, 'rx');
%             hold on;
%             subplot(2,1,2);
%             plot(i,vehicleStateCheck.v_x,'bo', i, vehicleState.v_x, 'rx');
%             hold on;
            
            % metdata update
            metadata.pathValidity(i) = 1;
            replanCounter = replanCounter + 1;
            if (replanCounter == 5)
                replan = 1;
                replanCounter = 0;
            end
        else
            scenario = [];
            % When there is a corruption in the original data, then set the
            % vehicle position to the closest position in the measurement
            [vehicleState, nearestIndex] = setVehiclePositionToNearest(segment_m, indexes, vehicleState);
            if (isempty(nearestIndex))
                disp('Measurement is corrupt, stopping simulation');
                break;
            end
            replan = 1;
            metadata.pathValidity(i) = 0;
            replanCounter = 0;
        end
        path(i,1) = vehicleState.X;
        path(i,2) = vehicleState.Y;
        intentionPath(i,1) = vehicleState.X - coefficients(1).coefficients(1)*sin(vehicleState.theta);
        intentionPath(i,2) = vehicleState.Y + coefficients(1).coefficients(1)*cos(vehicleState.theta);
        vehicleStateMemory{i} = vehicleState;
    end
end

function [c0, c1] = transformVideoData(c0Video,c1Video,origo, pose)
    d = ((origo(1)-pose(1))^2 + ...
                (origo(2)-pose(2))^2)^0.5;
    dtheta = pose(3) - origo(3);
    rho = [pose(1)-origo(1) pose(2)-origo(2)];
    alfa = atan(rho(2)/rho(1)); % slope of the difference vector
    % adjusting alfa based on planes (1st quarter plane needs no alignment
    % as angle starts from there. When rho(1) is zero,
    % angle needs no aligment, as rho(2) will determine +/- inf --> +/-
    % pi()/2. When rho(2) is zero, then alignment is needed with respect to
    % rho(1) sign, as atan(0) is always 0, however, we may need -pi().
    if (rho(1) < 0 && rho(2) > 0)
        % 2nd quarter-plane
        alfa = pi() + alfa; % rotating back from pi()
    elseif (rho(1) > 0 && rho(2) < 0)
        % fourth quarter-plane
        alfa = 2*pi() + alfa; % rotating back from 2pi()
    elseif (rho(1) < 0 && rho(2) < 0)
        alfa = pi() + alfa; % rotating on from pi()
    elseif (all(rho==0))
        alfa = 0; % vector is zero vector, hence no angle can be defined
    elseif (rho(2) == 0)
        if (rho(1)>0)
            alfa = 0;
        elseif (rho(1) < 0)
            alfa = -pi();
        end
    end
    displacementVector(1) = d * cos(alfa-origo(3));
    displacementVector(2) = d * sin(alfa-origo(3));
    c0 = c0Video - displacementVector(2);
    c1 = tan(atan(c1Video) - dtheta);
end

function [c, npLocal] = transformCoefficients(X,Y,theta, origo, pose)
    d = ((origo(1)-pose(1))^2 + ...
                (origo(2)-pose(2))^2)^0.5;
    dtheta = pose(3) - origo(3);
    rho = [pose(1)-origo(1) pose(2)-origo(2)];
    alfa = atan(rho(2)/rho(1)); % slope of the difference vector
    % adjusting alfa based on planes (1st quarter plane needs no alignment
    % as angle starts from there. When rho(1) is zero,
    % angle needs no aligment, as rho(2) will determine +/- inf --> +/-
    % pi()/2. When rho(2) is zero, then alignment is needed with respect to
    % rho(1) sign, as atan(0) is always 0, however, we may need -pi().
    if (rho(1) < 0 && rho(2) > 0)
        % 2nd quarter-plane
        alfa = pi() + alfa; % rotating back from pi()
    elseif (rho(1) > 0 && rho(2) < 0)
        % fourth quarter-plane
        alfa = 2*pi() + alfa; % rotating back from 2pi()
    elseif (rho(1) < 0 && rho(2) < 0)
        alfa = pi() + alfa; % rotating on from pi()
    elseif (all(rho==0))
        alfa = 0; % vector is zero vector, hence no angle can be defined
    elseif (rho(2) == 0)
        if (rho(1)>0)
            alfa = 0;
        elseif (rho(1) < 0)
            alfa = -pi();
        end
    end
    displacementVector(1) = d * cos(alfa-origo(3));
    displacementVector(2) = d * sin(alfa-origo(3));
    % doing the transformation
    T = [cos(dtheta) sin(dtheta); -sin(dtheta) cos(dtheta)];
    npLocal = ([X Y] - [displacementVector(1) displacementVector(2)])*T';
    thetaLocal = theta - dtheta;
    for i=1:size(npLocal,1)-1
        c(i).coefficients = OrderPolynomialRegression (npLocal(i,1), npLocal(i+1,1), npLocal(i,2), npLocal(i+1,2), thetaLocal(i), thetaLocal(i+1));
        c(i).sectionBorders = [npLocal(i,1), npLocal(i+1,1)];
    end
end

function vehicleState = initVehicleState(segment_m, indexes, index, mode)
if (mode=="kinematic")
    vehicleState.X = segment_m(index, indexes.X_abs);
    vehicleState.Y = segment_m(index, indexes.Y_abs);
    vehicleState.theta = segment_m(index, indexes.theta_calc);
    vehicleState.yawRate = segment_m(index, indexes.YawRate);
    vehicleState.vx = segment_m(index, indexes.VelocityX);
    vehicleState.vy = 0;
    vehicleState.ax = 0;
    vehicleState.ay = vehicleState.vx*vehicleState.yawRate;
    vehicleState.t0 = segment_m(index, indexes.Relative_time);    
    vehicleState.steeringAngle = 0;
    vehicleState.time = 0;
elseif (mode =="dynamic")
    lf = 1; lr = 1.5; r = 0.309725;
    
    vehicleState.X = segment_m(index, indexes.X_abs);
    vehicleState.Y = segment_m(index, indexes.Y_abs);
    vehicleState.theta = segment_m(index, indexes.theta_calc);
    vehicleState.v_x = segment_m(index, indexes.VelocityX);
    vehicleState.v_y = 0;
    vehicleState.t0 = segment_m(index, indexes.Relative_time);
    vehicleState.yawRate = segment_m(index, indexes.YawRate);
    vehicleState.steeringAngle = 0;
    vehicleState.v_fx = vehicleState.v_x;
    vehicleState.v_fy = vehicleState.v_y + vehicleState.yawRate*lf;
    vehicleState.v_rx = vehicleState.v_x;
    vehicleState.v_ry = vehicleState.v_y - vehicleState.yawRate*lr;
    vehicleState.drho_f = vehicleState.v_fx/r;
    vehicleState.drho_r = vehicleState.v_rx/r;
    vehicleState.v_fx_v = cos(vehicleState.steeringAngle)*vehicleState.v_fx + sin(vehicleState.steeringAngle)*vehicleState.v_fy;
    vehicleState.v_fy_v = -sin(vehicleState.steeringAngle)*vehicleState.v_fx + cos(vehicleState.steeringAngle)*vehicleState.v_fy;
    vehicleState.v_rx_v = vehicleState.v_rx;
    vehicleState.v_ry_v = vehicleState.v_ry;
    
    vehicleState.vx = (vehicleState.v_x^2 + vehicleState.v_y^2)^0.5;
    vehicleState.ay = vehicleState.v_x*vehicleState.yawRate;
    
    vehicleState.M_af = 0; vehicleState.M_ar = 0;
    vehicleState.M_bf = 0; vehicleState.M_br = 0;
    vehicleState.time = 0;
elseif (mode =="dynamicSimplified")
    lf = 1; lr = 1.5; r = 0.309725;
    
    vehicleState.X = segment_m(index, indexes.X_abs);
    vehicleState.Y = segment_m(index, indexes.Y_abs);
    vehicleState.theta = segment_m(index, indexes.theta_calc);
    vehicleState.v_x = segment_m(index, indexes.VelocityX);
    vehicleState.v_y = 0;
    vehicleState.t0 = segment_m(index, indexes.Relative_time);
    vehicleState.yawRate = segment_m(index, indexes.YawRate);
    vehicleState.steeringAngle = 0;
    vehicleState.v_fx = vehicleState.v_x;
    vehicleState.v_fy = vehicleState.v_y + vehicleState.yawRate*lf;
    vehicleState.v_rx = vehicleState.v_x;
    vehicleState.v_ry = vehicleState.v_y - vehicleState.yawRate*lr;
    vehicleState.drho_f = vehicleState.v_fx/r;
    vehicleState.drho_r = vehicleState.v_rx/r;
    vehicleState.v_fx_v = cos(vehicleState.steeringAngle)*vehicleState.v_fx + sin(vehicleState.steeringAngle)*vehicleState.v_fy;
    vehicleState.v_fy_v = -sin(vehicleState.steeringAngle)*vehicleState.v_fx + cos(vehicleState.steeringAngle)*vehicleState.v_fy;
    vehicleState.v_rx_v = vehicleState.v_rx;
    vehicleState.v_ry_v = vehicleState.v_ry;
    
    vehicleState.vx = (vehicleState.v_x^2 + vehicleState.v_y^2)^0.5;
    vehicleState.ay = vehicleState.v_x*vehicleState.yawRate;
    
    vehicleState.M_af = 0; vehicleState.M_ar = 0;
    vehicleState.M_bf = 0; vehicleState.M_br = 0;
end
end

function vehicleState = loadModel(vehicleState)
    b = 0; %0.0214;
    vehicleState.M_bf = vehicleState.v_fx * b;
    vehicleState.M_br = vehicleState.v_rx * b;
end

function [vehicleState, nearestIndex] = setVehiclePositionToNearest(segment_m, indexes, vehicleState)
global modelMode    
nearestIndex = getNearestIndex(segment_m(:,indexes.X_abs),segment_m(:,indexes.Y_abs),[vehicleState.X vehicleState.Y]);
    if (~isempty(nearestIndex))
        vehicleState = initVehicleState(segment_m, indexes, min(size(segment_m,1),nearestIndex(1,1)+1), modelMode);
    end
end

function scenario = cutScenario(segment_m, indexes, nearestIndex, vehicleState, window)
    % window = [pastDistance FutureDistance]
    vehiclePose = [vehicleState.X vehicleState.Y vehicleState.theta];
    distances = ((segment_m(:, indexes.X_abs) - vehiclePose(1)).^2+(segment_m(:, indexes.Y_abs) - vehiclePose(2)).^2).^0.5;
    if (window(1) == 0)
        scenarioStart = nearestIndex;
    else
        scenarioStart = find(distances(1:nearestIndex) <= window(1),1);
    end
    scenarioStop = find(distances(nearestIndex:end) > window(2),1);
    if (isempty(scenarioStart) || isempty(scenarioStop))
        scenario = [];
    else
        scenarioStop = nearestIndex + scenarioStop - 1;
        if (any(abs(diff(distances(scenarioStart:scenarioStop))) > 3))
            % there is discontinuity in the data
            scenario = [];
        else
            scenario = segment_m(scenarioStart:scenarioStop, :);
        end
    end
end

function [priorPath, u, dy] = planner (scenario, indexes, previousPosteriorPath, parameters, net)
    global modelID
    % some planner, which returns the priorPath
    % priorPath = NX2 array, path points ahead of the vehicle with a certain
    % step-size (1m)
    % priorPath = [scenario(:, indexes.X_abs) scenario(:, indexes.Y_abs)]; % temporary: ground-truth
    P_npDistances = parameters.P_npDistances;
    u = []; dy = [];
    switch modelID
        case "ldm"
            [priorPath, u, dy] = ldm (scenario, indexes, previousPosteriorPath, P_npDistances, parameters.P_LDM);
        case "eldm"
            [priorPath, u ,dy] = eldm (scenario, indexes, previousPosteriorPath, P_npDistances, parameters.P_ELDM);
        case "groundTruth"
            [priorPath, u, dy] = groundTruth (scenario, indexes, previousPosteriorPath,false);
        case "groundTruthFiltered"
            [priorPath, u, dy] = groundTruth (scenario, indexes, previousPosteriorPath, true);
        case "modifiedGroundTruth"
            [priorPath, u, dy] = modifiedGroundTruth (scenario, indexes, previousPosteriorPath);
        case "gp"
            [priorPath, u, dy] = gpModel (scenario, indexes, previousPosteriorPath);
        case "pp3"
            [priorPath, u, dy] = pp3Model (scenario, indexes, previousPosteriorPath);
        case "lrm"
            [priorPath, u, dy] = lrmModel (scenario, indexes, previousPosteriorPath);
        otherwise
            priorPath = [scenario(:,indexes.X_abs) scenario(:,indexes.Y_abs)];
    end
end

function [path, u, dy] = groundTruth (scenario, indexes, previousPath, filtered)
    global parameters
    % This function calculates the node points in front of the vehicle,
    % then applies the driver model to calculate the offset to the
    % mid-lane. Finally, a curve is fit onto the node points to result the
    % final path. The entire planning happens in the global UTF frame
    % Inputs:
    % [1] scenario - array with all signals, cut for the planning horizon 
    % [2] previous path - the previously planned trajectory, used for
    % initializing the new trajectory

    %% step 0: calculating the planner frame
    % finding the closest point in the scenario
    distances = ((scenario(:,indexes.X_abs)-previousPath(end,1)).^2+...
        (scenario(:,indexes.Y_abs) - previousPath(end,2)).^2).^0.5;
    nearestPoint = find(distances == min(distances),1);
    
    % estimation of orientation:
    if (size(previousPath,1) > 1)
        theta0 = atan2((previousPath(end,2)-previousPath(end-1,2)),(previousPath(end,1)-previousPath(end-1,1)));
    else
        theta0 = scenario(nearestPoint(1,1), indexes.theta_calc);
    end
    plannerFrame = [previousPath(end,:) theta0];
    Tplanner = [cos(theta0) sin(theta0); -sin(theta0) cos(theta0)];
    
    %% step 1: calculating the node points
    % Node Point Model
    corridorGlobalFrame = [scenario(nearestPoint(1,1):end,indexes.corrX) scenario(nearestPoint(1,1):end,indexes.corrY)];
    corridorPlannerFrame = (corridorGlobalFrame-plannerFrame(1:2))*Tplanner';
    [indeces, ~] = nodePointModel(corridorPlannerFrame, parameters.P_npDistances);
    
    if(filtered)
         referenceGlobalFrame = [scenario(nearestPoint(1,1):end,indexes.X_abs_mod) scenario(nearestPoint(1,1):end,indexes.Y_abs_mod)];
    else
         referenceGlobalFrame = [scenario(nearestPoint(1,1):end,indexes.X_abs) scenario(nearestPoint(1,1):end,indexes.Y_abs)];
    end
    referencePlannerFrame = (referenceGlobalFrame-plannerFrame(1:2))*Tplanner';
    Y_reference = referencePlannerFrame(indeces(2:end),2);
    
    % nominal road points
    X_nominal = corridorPlannerFrame(indeces(2:end),1);
    Y_nominal = corridorPlannerFrame(indeces(2:end),2);
    theta_nominal = scenario(indeces(2:end) + nearestPoint(1,1),indexes.LaneOrientation)+ ...
        scenario(indeces(2:end) + nearestPoint(1,1),indexes.theta_calc)- ...
        plannerFrame(3);
    
    %% step 2: calculating the offsets
    % Offset Model
    u = zeros(3,1);
    kappa_nominal = zeros(3,1);
    scenario(:,indexes.LaneCurvature) = movmean(scenario(:,indexes.LaneCurvature),100);
    for j=1:length(indeces)-1
        kappa_nominal(j,1) = mean(scenario(indeces(j) + nearestPoint(1,1):indeces(j+1) + nearestPoint(1,1),indexes.LaneCurvature));
    end
    u = kappa_nominal/0.001;
    
    % Reference calculation - transformed to local coordinate frame
    dy(1:length(indeces)-1,1) = (Y_reference - Y_nominal).*cos(theta_nominal);

    %% step 3: calculating final node points
    Y = Y_nominal + dy.*cos(theta_nominal);
    X = X_nominal - dy.*sin(theta_nominal);
    theta = theta_nominal;
    
    % extending with origin of the planner frame
    Y = [0; Y]; X = [0;X]; theta = [0; theta];
    
    %% step 4: curve fitting model
    refX = corridorPlannerFrame(:,1)';
    [pathPlannerFrame, ~] = trajectory_planner([X;Y;theta], indeces, refX',0);
    
    %% step 5: converting to global frame
    path = pathPlannerFrame*Tplanner + plannerFrame(1:2);
end

function [estimation, deviation, inputRaw] = gpGenerateEstimate(segment_m, indexes)
global parameters
    %% step 0: generate input data
    [~, ~, inputRaw, outputRaw] = prepareData(segment_m, indexes, parameters.PARAMS);
    
    %% step 1: parameter loading
    % read the learnt data from previous runs
    pathToParams = "C:\git\KDP\publications\GP\results\sparseGP_chosenParams\kpis";
    paramGP = dir(fullfile(pathToParams,"ETA_*"));

    for fileID = 1:length(paramGP)
        paramDriverID = str2num(paramGP(fileID).name(strfind(paramGP(fileID).name, 'driver_')+7:strfind(paramGP(fileID).name, '.')-1));
        if (paramDriverID == parameters.PARAMS.DriverID)
            paramData = load(fullfile(paramGP(fileID).folder, paramGP(fileID).name));
            for npID=1:10
                GP_params.hyp_opt_array{npID} = paramData.ETA(npID).hyp_opt;
                GP_params.output_estimation{npID} = paramData.ETA(npID).output_estimation;
                GP_params.input_estimation{npID} = paramData.ETA(npID).input_estimation;
            end
            c_in = paramData.ETA(npID).normFactors(1:7); % common for all GPs
            s_in = paramData.ETA(npID).normFactors(8:14); % common for all GPs
            c_out(1:10) = paramData.ETA(npID).normFactors(15); % one at each node point
            s_out = paramData.ETA(npID).normFactors(16:25); % one at each node point
            break;
        else
            paramData = [];
        end
    end

    %% step 3: norm and central
    for i=1:size(inputRaw,2)
        input(:,i) = (inputRaw(:,i)-c_in(i))/s_in(i);
    end
    for i=1:size(outputRaw,2)
        output(:,i) = (outputRaw(:,i)-c_out(i))/s_out(i);
    end

    %% step 4: generate output
    meanfunc = [];       % Start with a zero mean prior
    eval(strcat('covfunc = ',parameters.PARAMS.KERNEL_TYPE));    % Squared Exponental covariance function
    % ID problem            
    likfunc = @likGauss;    % Gaussian likelihood
    for i = 1:length(GP_params.hyp_opt_array)
        hyp = struct('mean', [], 'cov', 0, 'lik', -1);
        hyp.cov = GP_params.hyp_opt_array{i}.cov;
        hyp.lik = GP_params.hyp_opt_array{i}.lik;
        [estimationGP(:,i), deviationGP(:,i)] = gp(hyp, @infGaussLik, meanfunc, covfunc, likfunc, GP_params.input_estimation{i}, GP_params.output_estimation{i}, input); % extract the mean and covarance functions
    end

    %% step 5: node point calculation
    for i=1:size(estimationGP,2)
        estimation(:,i) = estimationGP(:,i)*s_out(i)+c_out(i);
        deviation(:,i) = deviationGP(:,i)*s_out(i)+c_out(i);
    end
end

function estimation = pp3GenerateEstimate(segment_m, indexes)
global parameters
    %% step 0: generate input data
    [~, ~, inputRaw, outputRaw] = prepareData(segment_m, indexes, parameters.PARAMS);
    
    %% step 1: parameter loading
    % read the learnt data from previous runs
    pathToParams = "C:\git\KDP\publications\GP\results\pp3data";
    paramGP = dir(fullfile(pathToParams,"KPI__driver*"));
    paramGP_ = dir(fullfile(pathToParams,"ETA__driver*"));
    for fileID = 1:length(paramGP)
        paramDriverID = str2num(paramGP(fileID).name(strfind(paramGP(fileID).name, 'driver_')+7:strfind(paramGP(fileID).name, '.')-1));
        if (paramDriverID == parameters.PARAMS.DriverID)
            paramData = load(fullfile(paramGP(fileID).folder, paramGP(fileID).name));
            for npID=1:10
                PP3_params.PP3_coeff{npID,1} = paramData.KPI{npID}(7:end-1);
            end
            break;
        else
            PP3_params = [];
        end
    end
    for fileID = 1:length(paramGP_)
        paramDriverID = str2num(paramGP_(fileID).name(strfind(paramGP_(fileID).name, 'driver_')+7:strfind(paramGP_(fileID).name, '.')-1));
        if (paramDriverID == parameters.PARAMS.DriverID)
            paramData = load(fullfile(paramGP_(fileID).folder, paramGP_(fileID).name));
            for npID=1:10
                PP3_params.output_estimation{npID} = paramData.ETA(npID).output_estimation;
                PP3_params.input_estimation{npID} = paramData.ETA(npID).input_estimation;
            end
            c_in = paramData.ETA(npID).normFactors(1:8); % common for all GPs
            s_in = paramData.ETA(npID).normFactors(9:16); % common for all GPs
            c_out(1:10) = paramData.ETA(npID).normFactors(17); % one at each node point
            s_out = paramData.ETA(npID).normFactors(18:27); % one at each node point
            break;
        end
    end

    %% step 3: norm and central
    for i=1:size(inputRaw,2)
        input(:,i) = (inputRaw(:,i)-c_in(i))/s_in(i);
    end
    for i=1:size(outputRaw,2)
        output(:,i) = (outputRaw(:,i)-c_out(i))/s_out(i);
    end

    %% step 4: generate output
    for i = 1:length(PP3_params.PP3_coeff)
        estimationPP3(:,i) = reFitPolynomial(input, PP3_params.PP3_coeff{npID}', 3);
    end

    %% step 5: node point calculation
    for i=1:size(estimationPP3,2)
        estimation(:,i) = estimationPP3(:,i)*s_out(i)+c_out(i);
    end
end

function y_ = reFitPolynomial(X, c, p)
y_ = zeros(size(X,1),1);
for i=1:p
    y_ = y_+(c((i-1)*size(X,2)+1:i*size(X,2))'*(X.^i)')';
end
end

function [path, u, dy] = gpModel (scenario, indexes, previousPath)
    global parameters
    % This function calculates the node points in front of the vehicle,
    % then applies the driver model to calculate the offset to the
    % mid-lane. Finally, a curve is fit onto the node points to result the
    % final path. The entire planning happens in the global UTF frame
    % Inputs:
    % [1] scenario - array with all signals, cut for the planning horizon 
    % [2] previous path - the previously planned trajectory, used for
    % initializing the new trajectory
    u = []; dy=[];

    %% step 1: path planning
    % finding the closest point in the scenario
    distances = ((scenario(:,indexes.X_abs)-previousPath(end,1)).^2+...
        (scenario(:,indexes.Y_abs) - previousPath(end,2)).^2).^0.5;
    nearestPoint = find(distances == min(distances),1);
    
    % estimation of orientation:
    if (size(previousPath,1) > 1)
        theta0 = atan2((previousPath(end,2)-previousPath(end-1,2)),(previousPath(end,1)-previousPath(end-1,1)));
    else
        theta0 = scenario(nearestPoint(1,1), indexes.theta_calc);
    end
    plannerFrame = [previousPath(end,:) theta0];
    Tplanner = [cos(theta0) sin(theta0); -sin(theta0) cos(theta0)];
    
    %% step 2: calculating the node points
    % Node Point Model
    corridorGlobalFrame = [scenario(nearestPoint(1,1):end,indexes.corrX) scenario(nearestPoint(1,1):end,indexes.corrY)];
    corridorPlannerFrame = (corridorGlobalFrame-plannerFrame(1:2))*Tplanner';
    [indeces, ~] = nodePointModel(corridorPlannerFrame, [parameters.PARAMS.OUTPUT_SHIFT(1) diff(parameters.PARAMS.OUTPUT_SHIFT)]/250);
    
    referenceGlobalFrame = [scenario(nearestPoint(1,1):end,indexes.X_abs) scenario(nearestPoint(1,1):end,indexes.Y_abs)];
    referencePlannerFrame = (referenceGlobalFrame-plannerFrame(1:2))*Tplanner';
    Y_reference = referencePlannerFrame(indeces(2:end),2);
    
    % nominal road points
    X_nominal = corridorPlannerFrame(indeces(2:end),1);
    Y_nominal = corridorPlannerFrame(indeces(2:end),2);
    theta_nominal = scenario(indeces(2:end) + nearestPoint(1,1),indexes.LaneOrientation)+ ...
        scenario(indeces(2:end) + nearestPoint(1,1),indexes.theta_calc)- ...
        plannerFrame(3);
    
    % calculating final node points
    Y = Y_nominal + scenario(1,indexes.GP_1:indexes.GP_10)'.*cos(theta_nominal);
    X = X_nominal - scenario(1,indexes.GP_1:indexes.GP_10)'.*sin(theta_nominal);
    theta = theta_nominal;
    
    % extending with origin of the planner frame
    Y = [0; Y]; X = [0;X]; theta = [0; theta];
    
    %% step 4: curve fitting model
    refX = corridorPlannerFrame(:,1)';
    pathY = zeros(size(refX));
    n = 3;
    n = min(n, length(X)-1);
    c = polyfit(X(2:end),Y(2:end),n);
    for i=1:n+1
        pathY = pathY + c(i)*refX.^(n-(i-1));
    end
    pathPlannerFrame = [refX' pathY'];
    %[pathPlannerFrame, ~] = trajectory_planner([X;Y;theta], indeces, refX',0);
    
    %% step 5: converting to global frame
    path = pathPlannerFrame*Tplanner + plannerFrame(1:2);
end

function [path, u, dy] = pp3Model (scenario, indexes, previousPath)
    global parameters
    % This function calculates the node points in front of the vehicle,
    % then applies the driver model to calculate the offset to the
    % mid-lane. Finally, a curve is fit onto the node points to result the
    % final path. The entire planning happens in the global UTF frame
    % Inputs:
    % [1] scenario - array with all signals, cut for the planning horizon 
    % [2] previous path - the previously planned trajectory, used for
    % initializing the new trajectory
    u = []; dy=[];

    %% step 1: path planning
    % finding the closest point in the scenario
    distances = ((scenario(:,indexes.X_abs)-previousPath(end,1)).^2+...
        (scenario(:,indexes.Y_abs) - previousPath(end,2)).^2).^0.5;
    nearestPoint = find(distances == min(distances),1);
    
    % estimation of orientation:
    if (size(previousPath,1) > 1)
        theta0 = atan2((previousPath(end,2)-previousPath(end-1,2)),(previousPath(end,1)-previousPath(end-1,1)));
    else
        theta0 = scenario(nearestPoint(1,1), indexes.theta_calc);
    end
    plannerFrame = [previousPath(end,:) theta0];
    Tplanner = [cos(theta0) sin(theta0); -sin(theta0) cos(theta0)];
    
    %% step 2: calculating the node points
    % Node Point Model
    corridorGlobalFrame = [scenario(nearestPoint(1,1):end,indexes.corrX) scenario(nearestPoint(1,1):end,indexes.corrY)];
    corridorPlannerFrame = (corridorGlobalFrame-plannerFrame(1:2))*Tplanner';
    [indeces, ~] = nodePointModel(corridorPlannerFrame, [parameters.PARAMS.OUTPUT_SHIFT(1) diff(parameters.PARAMS.OUTPUT_SHIFT)]/250);
    
    referenceGlobalFrame = [scenario(nearestPoint(1,1):end,indexes.X_abs) scenario(nearestPoint(1,1):end,indexes.Y_abs)];
    referencePlannerFrame = (referenceGlobalFrame-plannerFrame(1:2))*Tplanner';
    Y_reference = referencePlannerFrame(indeces(2:end),2);
    
    % nominal road points
    X_nominal = corridorPlannerFrame(indeces(2:end),1);
    Y_nominal = corridorPlannerFrame(indeces(2:end),2);
    theta_nominal = scenario(indeces(2:end) + nearestPoint(1,1),indexes.LaneOrientation)+ ...
        scenario(indeces(2:end) + nearestPoint(1,1),indexes.theta_calc)- ...
        plannerFrame(3);
    
    % calculating final node points
    Y = Y_nominal + scenario(1,indexes.PP3_1:indexes.PP3_10)'.*cos(theta_nominal);
    X = X_nominal - scenario(1,indexes.PP3_1:indexes.PP3_10)'.*sin(theta_nominal);
    theta = theta_nominal;
    
    % extending with origin of the planner frame
    Y = [0; Y]; X = [0;X]; theta = [0; theta];
    
    %% step 4: curve fitting model
    refX = corridorPlannerFrame(:,1)';
    [pathPlannerFrame, ~] = trajectory_planner([X;Y;theta], indeces, refX',0);
    
    %% step 5: converting to global frame
    path = pathPlannerFrame*Tplanner + plannerFrame(1:2);
end

function [path, u, dy] = lrmModel (scenario, indexes, previousPath)
    global parameters
    % This function calculates the node points in front of the vehicle,
    % then applies the driver model to calculate the offset to the
    % mid-lane. Finally, a curve is fit onto the node points to result the
    % final path. The entire planning happens in the global UTF frame
    % Inputs:
    % [1] scenario - array with all signals, cut for the planning horizon 
    % [2] previous path - the previously planned trajectory, used for
    % initializing the new trajectory
    u = []; dy=[];

    %% step 1: path planning
    % finding the closest point in the scenario
    distances = ((scenario(:,indexes.X_abs)-previousPath(end,1)).^2+...
        (scenario(:,indexes.Y_abs) - previousPath(end,2)).^2).^0.5;
    nearestPoint = find(distances == min(distances),1);
    
    % estimation of orientation:
    if (size(previousPath,1) > 1)
        theta0 = atan2((previousPath(end,2)-previousPath(end-1,2)),(previousPath(end,1)-previousPath(end-1,1)));
    else
        theta0 = scenario(nearestPoint(1,1), indexes.theta_calc);
    end
    plannerFrame = [previousPath(end,:) theta0];
    Tplanner = [cos(theta0) sin(theta0); -sin(theta0) cos(theta0)];
    
    %% step 2: calculating the node points
    % Node Point Model
    corridorGlobalFrame = [scenario(nearestPoint(1,1):end,indexes.corrX) scenario(nearestPoint(1,1):end,indexes.corrY)];
    corridorPlannerFrame = (corridorGlobalFrame-plannerFrame(1:2))*Tplanner';
    [indeces, ~] = nodePointModel(corridorPlannerFrame, [parameters.PARAMS.OUTPUT_SHIFT(1) diff(parameters.PARAMS.OUTPUT_SHIFT)]/250);
    
    referenceGlobalFrame = [scenario(nearestPoint(1,1):end,indexes.X_abs) scenario(nearestPoint(1,1):end,indexes.Y_abs)];
    referencePlannerFrame = (referenceGlobalFrame-plannerFrame(1:2))*Tplanner';
    Y_reference = referencePlannerFrame(indeces(2:end),2);
    
    % nominal road points
    X_nominal = corridorPlannerFrame(indeces(2:end),1);
    Y_nominal = corridorPlannerFrame(indeces(2:end),2);
    theta_nominal = scenario(indeces(2:end) + nearestPoint(1,1),indexes.LaneOrientation)+ ...
        scenario(indeces(2:end) + nearestPoint(1,1),indexes.theta_calc)- ...
        plannerFrame(3);
    
    % calculating final node points
    Y = Y_nominal + scenario(1,indexes.LRM_1:indexes.LRM_10)'.*cos(theta_nominal);
    X = X_nominal - scenario(1,indexes.LRM_1:indexes.LRM_10)'.*sin(theta_nominal);
    theta = theta_nominal;
    
    % extending with origin of the planner frame
    Y = [0; Y]; X = [0;X]; theta = [0; theta];
    
    %% step 4: curve fitting model
    refX = corridorPlannerFrame(:,1)';
    [pathPlannerFrame, ~] = trajectory_planner([X;Y;theta], indeces, refX',0);
    
    %% step 5: converting to global frame
    path = pathPlannerFrame*Tplanner + plannerFrame(1:2);
end

function [dataOut, c,s]= normAndCentral(dataIn)
    for i=1:size(dataIn,2)
        c(i) = mean(dataIn(:,i));
        s(i) = std(dataIn(:,i));
        dataOut(:,i) = (dataIn(:,i)-mean(dataIn(:,i)))/std(dataIn(:,i));
    end
end

function [path, u, dy] = modifiedGroundTruth (scenario, indexes, previousPath)
    global parameters
    % This function calculates the node points in front of the vehicle,
    % then applies the driver model to calculate the offset to the
    % mid-lane. Finally, a curve is fit onto the node points to result the
    % final path. The entire planning happens in the global UTF frame
    % Inputs:
    % [1] scenario - array with all signals, cut for the planning horizon 
    % [2] previous path - the previously planned trajectory, used for
    % initializing the new trajectory

    %% step 0: calculating the planner frame
    % finding the closest point in the scenario
    distances = ((scenario(:,indexes.X_abs_mod)-previousPath(end,1)).^2+...
        (scenario(:,indexes.Y_abs_mod) - previousPath(end,2)).^2).^0.5;
    nearestPoint = find(distances == min(distances),1);
    
    % estimation of orientation:
    if (size(previousPath,1) > 1)
        theta0 = atan2((previousPath(end,2)-previousPath(end-1,2)),(previousPath(end,1)-previousPath(end-1,1)));
    else
        theta0 = scenario(nearestPoint(1,1), indexes.theta_calc);
    end
    plannerFrame = [previousPath(end,:) theta0];
    Tplanner = [cos(theta0) sin(theta0); -sin(theta0) cos(theta0)];
    
    %% step 1: calculating the node points
    % Node Point Model
    corridorGlobalFrame = [scenario(nearestPoint(1,1):end,indexes.corrX) scenario(nearestPoint(1,1):end,indexes.corrY)];
    corridorPlannerFrame = (corridorGlobalFrame-plannerFrame(1:2))*Tplanner';
    [indeces, ~] = nodePointModel(corridorPlannerFrame, parameters.P_npDistances);
    
    referenceGlobalFrame = [scenario(nearestPoint(1,1):end,indexes.X_abs_mod) scenario(nearestPoint(1,1):end,indexes.Y_abs_mod)];
    referencePlannerFrame = (referenceGlobalFrame-plannerFrame(1:2))*Tplanner';
    Y_reference = referencePlannerFrame(indeces(2:end),2);
    
    % nominal road points
    X_nominal = corridorPlannerFrame(indeces(2:end),1);
    Y_nominal = corridorPlannerFrame(indeces(2:end),2);
    theta_nominal = scenario(indeces(2:end) + nearestPoint(1,1),indexes.LaneOrientation)+ ...
        scenario(indeces(2:end) + nearestPoint(1,1),indexes.theta_calc)- ...
        plannerFrame(3);
    
    %% step 2: calculating the offsets
    % Offset Model
    u = zeros(3,1);
    kappa_nominal = zeros(3,1);
    scenario(:,indexes.LaneCurvature) = movmean(scenario(:,indexes.LaneCurvature),100);
    for j=1:length(indeces)-1
        kappa_nominal(j,1) = mean(scenario(indeces(j) + nearestPoint(1,1):indeces(j+1) + nearestPoint(1,1),indexes.LaneCurvature));
    end
    % scaling of the predictor variables
    kappa_lim = 0.001;
%     u_lim = [ones(6,1)*kappa_lim; 1];
%     u = [max(kappa_nominal,0); min(kappa_nominal,0); 1];
    u_lim = ones(3,1)*kappa_lim;
    u = kappa_nominal;
    u = u./u_lim;
    
    % Reference calculation - transformed to local coordinate frame
    dy(1:3,1) = (Y_reference - Y_nominal).*cos(theta_nominal);

    %% step 3: calculating final node points
    Y = Y_nominal + dy.*cos(theta_nominal);
    X = X_nominal - dy.*sin(theta_nominal);
    theta = theta_nominal;
    
    % extending with origin of the planner frame
    Y = [0; Y]; X = [0;X]; theta = [0; theta];
    
    %% step 4: curve fitting model
    refX = corridorPlannerFrame(:,1)';
    [pathPlannerFrame, ~] = trajectory_planner([X;Y;theta], indeces, refX',0);
    
    %% step 5: converting to global frame
    path = pathPlannerFrame*Tplanner + plannerFrame(1:2);
end

function [targetSteeringAngle] = controller(path, vehicleState)
% pure-pursuit
p_lookAheadTime = 0.34;
p_wheelBase = 2.7; % meter

lad = p_lookAheadTime*vehicleState.vx;
% convert path to local frame
T = [cos(vehicleState.theta) sin(vehicleState.theta); ...
    -sin(vehicleState.theta) cos(vehicleState.theta)];
localPath = (path - [vehicleState.X vehicleState.Y])*T';
localPath(localPath(:,1)<0,2) = 0;
localPath(localPath(:,1)<0,1) = 0;
% find look ahead point
idx = find((localPath(:,1).^2+localPath(:,2).^2).^0.5 >= lad,1);
L = localPath(idx(1,1),1);
y = spline(localPath(localPath(:,1)>0,1), localPath(localPath(:,1)>0,2), L);

% pure-pursuit relation
targetSteeringAngle = atan(p_wheelBase*(2*y)/L^2);
end

function [targetSteeringAngle] = controllerLite(coefficients, vehicleState, dT)
% pure-pursuit
p_lookAheadTime = 0.35;
p_wheelBase = 2.7; % meter

lad = p_lookAheadTime*vehicleState.vx;
% selecting right polyLine
c = coefficients(end).coefficients;
for i=1:length(coefficients)
    if (coefficients(i).sectionBorders(1) <= lad && coefficients(i).sectionBorders(2)>lad)
        c = coefficients(i).coefficients;
        break;
    end
end
y = c(1)+c(2)*lad+c(3)*lad^2+c(4)*lad^3;
L = (lad^2+y^2)^0.5;

% pure-pursuit relation
targetSteeringAngle = atan(p_wheelBase*(2*y)/L^2);


%% MPC mode
global x_k1
% prepare input
c1 = coefficients(1);
c2 = coefficients(2);
c3 = coefficients(3);

% resconstructing path local
dx = vehicleState.vx*dT;
x1 = 0:dx:c1.sectionBorders(2)-dx;
y1 = c1.coefficients(1) + c1.coefficients(2)*x1 +  c1.coefficients(3)*x1.^2 +  c1.coefficients(4)*x1.^3;
theta1 = atan(c1.coefficients(2) +  2*c1.coefficients(3)*x1 +  3*c1.coefficients(4)*x1.^2);
x2 = c1.sectionBorders(2):dx:c2.sectionBorders(2);
y2 = c2.coefficients(1) + c2.coefficients(2)*x2 +  c2.coefficients(3)*x2.^2 +  c2.coefficients(4)*x2.^3;
theta2 = atan(c2.coefficients(2) +  2*c2.coefficients(3)*x2 +  3*c2.coefficients(4)*x2.^2);
x3 = c2.sectionBorders(2):dx:c3.sectionBorders(2);
y3 = c3.coefficients(1) + c3.coefficients(2)*x3 +  c3.coefficients(3)*x3.^2 +  c3.coefficients(4)*x3.^3;
theta3 = atan(c3.coefficients(2) +  2*c3.coefficients(3)*x3 +  3*c3.coefficients(4)*x3.^2);
pathLocal = [[x1 x2 x3]' [y1 y2 y3]'];
pathOrientation = [theta1 theta2 theta3]';

load("C:\git\KDP_Igneczi\publikációk\LaneWandering\data\MOE.mat");
if (pathLocal(1,2) <= -MOE(1,1))
    % left side drift reached intervention point
    q = 0.0037;
    rw = 47.62;
elseif (pathLocal(1,2)>= -MOE(1,2))
    q = 0.0031;
    rw = 44.92;
elseif (pathLocal(1,2) < 0)
    % left side drift
    q = 0.0043;
    rw = 72.38;
else
    q = 0.0055;
    rw = 80.28;
end

x = [0; 0; vehicleState.vx; 0; 0];
xa = [x-x_k1; 0; 0];

targetAcceleration = mpc(pathLocal, pathOrientation, 90, 90, rw, [q 0], 2.7, dT,  x_k1, xa, vehicleState.ay, vehicleState.steeringAngle, 0.14, 10, vehicleState.ay);
targetSteeringAngle = atan(targetAcceleration/vehicleState.vx^2 * 2.7);

x_k1 = x;
end

function u = mpc(pathLocal, pathOrientation, Np, Nc, rw, q, L, Ts,  x_k1, xa, aeta0, delta0, deltamax, ay_max, u_k1)
%% INTRODUCTION
% This function is the core MPC algorithm part.
% created by Gergo Igneczi @ Vehicle Research Center of Szechenyi Istvan
% University

if (x_k1(3) < 3)
    % speed is low, MPC will not work
    u = 0;
else
    dim = 2; % number of outputs
    % producing prediction matrices
    pred_matrix = zeros(min(dim*Np,1000),1);

    x_points = pathLocal(:,1);
    y_points = pathLocal(:,2);

    for i = 1:Np
        %pred_matrix(dim*i-(dim-1),1) = x_points(i);
        pred_matrix(dim*i-(dim-1),1) = y_points(i);
        pred_matrix(dim*i-(dim-2),1) = pathOrientation(i);
    end
    Rs_rk = pred_matrix;

    I = eye(min(Nc,1000));
    %% determining state constraints    

    u_max = ay_max;
    u_min =  ay_max;

    %% initializing helper matrices
    M = zeros(min(1000,2*size(u_k1,1)*Nc),min(2*Nc,1000));
    for (i=1:2*size(u_k1,1)*Nc)
        if (i<=(Nc))
            k=1;
            while ((2*k-1)/2 <= i)
                M (i,2*k-1) = 1;
                k = k + 1;
            end
        elseif (i <=(Nc*2))
            k=1;
            while ((2*k-1)/2 <= i-Nc)
                M (i,2*k-1) = -1;
                k = k + 1;
            end
        elseif (i<=(Nc*3))
            k=1;
            while (2*k/2 <= i-size(x_k1,1)*Nc/2)
                M (i,2*k) = 1;
                k = k+1;
            end
        else
            k=1;
            while (2*k/2 <= i-size(x_k1,1)*Nc/2-Nc)
                M (i,2*k) = -1;
                k = k+1;
            end
        end
    end


    %% updating state matrices
    %Ad = [1 0 Ts 0 0; 0 1 0 Ts 0; 0 0 1 0 -Ts/L*x_k1(3)^2*tan(delta0); 0 0 0 1 0; 0 0 0 0 1];
    %Bd = [0 0 0 Ts/L*x_k1(3)^2*1/(cos(delta0))^2 Ts/L*x_k1(3)*1/(cos(delta0))^2]';
    Ad = [1 0 Ts 0 0; 0 1 0 Ts 0; 0 0 1 0 -aeta0*Ts; 0 0 0 1 0;0 0 0 0 1];
    Bd = [0 0 0 Ts Ts/x_k1(3)]';
    Cd = [0 1 0 0 0; 0 0 0 0 1];
    
    n = size(Ad,1); m = size(Cd,1); k = 1; %size(Bd,2);
    
    %% augmented model
    A = [Ad zeros(n,m); Cd*Ad eye(m,m)];
    B = [Bd; Cd*Bd];
    C = [zeros(m,n) eye(m,m)];
    %% matrix generation
    F = zeros(Np*m,m+n);
    for i=1:Np
        if (m>1)
            F(m*i-(m-1):m*i,1:(m+n))=C*A^i;
        else
            F(i,1:(m+n))=C*A^i;
        end
    end
    S = zeros(m*Np,Nc);
    for i=1:Np
        for j=1:Nc
            if(j>i)
                S(m*i-(m-1):m*i,j)=zeros(m,k);
            else
                S(m*i-(m-1):m*i,j)=C*A^((i-1)-(j-1))*B;
            end
        end
    end

    %% control calculation
    % Unconstrained results
    dU = zeros(min(1000,k*Nc),1);
    %dU = inv(S'*S+rw*I)*(S'*Rs_rk-S'*F*xa);
    R = rw*I;
    Q = zeros(m*Np, m*Np);
    for i=1:Np
        Q(i*numel(q)-(numel(q)-1):i*numel(q),i*numel(q)-(numel(q)-1):i*numel(q)) = diag(q);
    end
    dU = inv(S'*Q*S+R)*S'*Q*(Rs_rk-F*xa);
 
    % Constrained results
    gamma = zeros(min(1000,2*size(u_k1,1)*Nc),1);

    for i=1:2*size(u_k1,1)*Nc
        if (i<=Nc)
            gamma(i,1) = u_max(1)-u_k1(1);
        elseif (i<=2*Nc)
            gamma(i,1) = u_min(1)+u_k1(1);
        end
    end

    if (all(M(:,1:2:end)*dU<=gamma))
        %do nothing
    else
        %Solving Hildreth's QP problem
        E = (S'*S+rw*I)*2;
        F_ = -2*S'*(Rs_rk-F*xa);
        H = E;
        f = F_;
        A_cons = M(:,1:2:end);
        b = gamma;
        eta = x_k1;
        [n1,m1]=size(A_cons);
        eta=-H\f;
        kk=0;
        for i=1:n1
            if (A_cons(i,:)*eta>b(i)) 
                kk=kk+1;
            else
                kk=kk+0;
            end
        end
        if (kk==0) 
            % do nothing 
        else
            P=A_cons*(H\A_cons');
            d=(A_cons*(H\f)+b);
            [n,m]=size(d);
            x_ini=zeros(n,m);
            lambda=x_ini;
            al=10;
            for km=1:38
                %find the elements in the solution vector one by one
                % km could be larger if the Lagranger multiplier has a slow
                % convergence rate.
                lambda_p=lambda;
                for i=1:n
                    w= P(i,:)*lambda-P(i,i)*lambda(i,1);
                    w=w+d(i,1);
                    la=-w/P(i,i);
                    lambda(i,:)=max(0,la);
                end
                al=(lambda-lambda_p)'*(lambda-lambda_p);
                if (al<10e-8)
                    break; 
                end
            end
            dU=-H\f -H\A_cons'*lambda;
        end
    end

    u = u_k1+dU(1);
    Y = F*xa+S*dU;
    %u(1) = min(max(u(1),-ay_max),ay_max);
%     plot(x_points,y_points,x_points(1:Np),Y(1:2:end), 'LineWidth', 2, 'LineStyle', '--');
%     grid on; ylim([-15, 15]);
%     pause(0.1);

end


end

function vehicleState = speedController(vehicleState)
    setSpeed = 35 / 3.6;
    P = 0.1;
    speedError = setSpeed - vehicleState.vx;
    vehicleState.M_av = 0; %speedError * P;
end

function scenarioFinished = scenarioFinishChecker(path, vehicleState)
    p_lookAheadTime = 1;

    lad = p_lookAheadTime*vehicleState.vx;
    % convert path to local frame
    T = [cos(vehicleState.theta) sin(vehicleState.theta); ...
    -sin(vehicleState.theta) cos(vehicleState.theta)];
    localPath = (path - [vehicleState.X vehicleState.Y])*T';
    localPath(localPath(:,1)<0,2) = 0;
    localPath(localPath(:,1)<0,1) = 0;
    % find look ahead point
    idx = find((localPath(:,1).^2+localPath(:,2).^2).^0.5 >= lad,1);
    if (isempty(idx))
        scenarioFinished = 1;
    else
        scenarioFinished = 0;
    end
end

function [path, U, dy] = ldm (scenario, indexes, previousPath, P_npDistances, P_LDM)
global dt
    % This function calculates the node points in front of the vehicle,
    % then applies the driver model to calculate the offset to the
    % mid-lane. Finally, a curve is fit onto the node points to result the
    % final path. The entire planning happens in the global UTF frame
    % Inputs:
    % [1] scenario - array with all signals, cut for the planning horizon 
    % [2] previous path - the previously planned trajectory, used for
    % initializing the new trajectory
    % Parameters:
    %% step 0: calculating the planner frame
    % finding the closest point in the scenario
    tic;
    N = length(P_npDistances); % number of node points

    distances = ((scenario(:,indexes.X_abs)-previousPath(end,1)).^2+...
        (scenario(:,indexes.Y_abs) - previousPath(end,2)).^2).^0.5;
    nearestPoint = find(distances == min(distances),1);
    
    scenario(:,indexes.LaneCurvature) = movmean(scenario(:,indexes.LaneCurvature),100);
    
    % estimation of orientation:
    if (size(previousPath,1) > 1)
        theta0 = atan2((previousPath(end,2)-previousPath(end-1,2)),(previousPath(end,1)-previousPath(end-1,1)));
    else
        theta0 = scenario(nearestPoint(1,1), indexes.theta_calc);
    end
    plannerFrame = [previousPath(end,:) theta0];
    Tplanner = [cos(theta0) sin(theta0); -sin(theta0) cos(theta0)];
    
    %% step 1: calculating the node points
    % Node Point Model
    corridorGlobalFrame = [scenario(nearestPoint(1,1):end,indexes.corrX) scenario(nearestPoint(1,1):end,indexes.corrY)];
    corridorPlannerFrame = (corridorGlobalFrame-plannerFrame(1:2))*Tplanner';
    [indeces, ~] = nodePointModel(corridorPlannerFrame, P_npDistances);
    
    % nominal road points
    X_nominal = corridorPlannerFrame(indeces(2:end),1);
    Y_nominal = corridorPlannerFrame(indeces(2:end),2);
    theta_nominal = scenario(indeces(2:end) + nearestPoint(1,1),indexes.LaneOrientation)+ ...
        scenario(indeces(2:end) + nearestPoint(1,1),indexes.theta_calc)- ...
        plannerFrame(3);
    
    %% step 2: calculating the offsets
    % Offset Model
    U = zeros(N,1);
    for j=1:length(indeces)-1
        U(j,1) = mean(scenario(indeces(j) + nearestPoint(1,1):indeces(j+1) + nearestPoint(1,1),indexes.LaneCurvature));
    end
    % scaling of the predictor variables
    kappa_lim = 0.001;
    U_lim = ones(N,1)*kappa_lim;
    U = U./U_lim;

    % The ouput of the model is the linear combination of predictor variables
    dy = P_LDM*U;
    dy(1:N) = min(max(-1.25,dy(1:N)), 1.25);
    
    %% step 3: calculating final node points
    Y = Y_nominal + dy.*cos(theta_nominal);
    X = X_nominal - dy.*sin(theta_nominal);
    theta = theta_nominal;
    
    % extending with origin of the planner frame
    Y = [0; Y]; X = [0;X]; theta = [0; theta];
    
    %% step 4: curve fitting model
    refX = corridorPlannerFrame(:,1)';
    [pathPlannerFrame, ~] = trajectory_planner([X;Y;theta], indeces, refX',0);
    
    %% step 5: converting to global frame
    path = pathPlannerFrame*Tplanner + plannerFrame(1:2);
    dt = toc;
end

function [c_3polynomials, U, dy, X,Y,theta] = eldmLite (scenario, indexes, vehicleState, coefficients_k1, P_npDistances, P_ELDM, c0, c1)
    % This function calculates the node points in front of the vehicle,
    % then applies the driver model to calculate the offset to the
    % mid-lane. Finally, a curve is fit onto the node points to result the
    % final path. The entire planning happens in the global UTF frame
    % Inputs:
    % [1] scenario - array with all signals, cut for the planning horizon 
    % [2] previous path - the previously planned trajectory, used for
    % initializing the new trajectory
    % Parameters:
    N = parameters.numberOfNodePoints; % number of node points
    %% step 0: Model input reconstruction
    % building U vector with average kappa values
    % scenario starts at ego position - NOTE: the commented solution is
    % applicable when the corridor is represented by a polynomial with at
    % least third order
    initLastTraj = true;
    u = zeros(N,1);
%     for i=1:N
%         % c2 is scaled to curvature
%         if (i==1)
%             u(i,1) = 2*scenario(1,indexes.LaneCurvature)+3*scenario(1,indexes.c3)*P_npDistances(i);
%         else
%             u(i,1) = 2*scenario(1,indexes.LaneCurvature)+3*scenario(1,indexes.c3)*(P_npDistances(i)+P_npDistances(i-1));
%         end            
%     end
    % Solution B
    % transforming corridor to the ego frame
    T = [cos(vehicleState.theta) sin(vehicleState.theta);
        -sin(vehicleState.theta) cos(vehicleState.theta)];
    corridorGlobalFrame = [scenario(:,indexes.corrX) scenario(:,indexes.corrY)];
    corridorEgoFrame = (corridorGlobalFrame - [vehicleState.X vehicleState.Y])*T';
    [indeces, ~] = nodePointModel(corridorEgoFrame, P_npDistances);
    for j=1:length(indeces)-1
        u(j,1) = mean(scenario(indeces(j):indeces(j+1),indexes.LaneCurvature));
    end
    
    % transforming due to E-LDM structure
    U = [max(u,0); min(u,0); 1];
    
    %% step 1: node point model
    % producing nominal values
    P_npDistances = P_npDistances*250;
    P_npDistances(1) = max(10,P_npDistances(1));
    for i=2:length(P_npDistances)
        P_npDistances(i) = max(10,P_npDistances(i-1)+P_npDistances(i));
    end

    nominalSelect = "";
    
    if (nominalSelect == "groundTruth")
        X_nominal = corridorEgoFrame(indeces(2:end),1);
        Y_nominal = corridorEgoFrame(indeces(2:end),2);
        theta_nominal = atan(scenario(indeces(2:end),indexes.LaneOrientation))+scenario(indeces(2:end),indexes.theta_calc)-vehicleState.theta;
    else
        for i=1:N
            X_nominal(i) = P_npDistances(i);
            Y_nominal(i) = c0 + c1*P_npDistances(i) + ...
                scenario(1,indexes.LaneCurvature)/2*P_npDistances(i)^2+ ...
                scenario(1,indexes.c3)/6*P_npDistances(i)^3;
            theta_nominal(i) = atan(c1 + ...
                scenario(1,indexes.LaneCurvature)*P_npDistances(i)+ ...
                scenario(1,indexes.c3)/2*P_npDistances(i)^2);
        end
    end
            
    
    %% step 2: offset model   
    % The ouput of the model is the linear combination of predictor variables
    dy = P_ELDM*U;
    % Limitation of the offset
    dy(1:3) = min(max(-1.25,dy(1:3)), 1.25);
    
    %% step 3: constructing node points
    % N+1 node points: 1 as the last point of the previous trajectory and N
    % from nominal points + offset
    X = zeros(N+1,1); Y = zeros(N+1,1); theta = zeros(N+1,1);
    % now find the valid polyline
    if (~isempty(coefficients_k1) && initLastTraj == true)
        c = coefficients_k1(N-1).coefficients;
        for i=1:N-1
            if (coefficients_k1(i).sectionBorders(1) <=0 && coefficients_k1(i).sectionBorders(2) > 0)
                c = coefficients_k1(i).coefficients;
                break;
            end
        end
    else
        c = zeros(1,4);
    end
    X(1) = 0; Y(1) = c(1); theta(1) = atan(c(2));
    for i=1:N
        X(i+1) = X_nominal(i) - sin(theta_nominal(i))*dy(i); % offset is always perpendicular to the refLine
        Y(i+1) = Y_nominal(i) + cos(theta_nominal(i))*dy(i);
        theta(i+1) = theta_nominal(i);
    end
    
    %% step 4: curve fitting
    for i=1:N
        c_3polynomials(i).coefficients = OrderPolynomialRegression (X(i), X(i+1), Y(i), Y(i+1), theta(i), theta(i+1));
        c_3polynomials(i).sectionBorders = [X(i), X(i+1)];
    end
end

function c = fit3rdorderPolynomial(X,Y)
    A = [ones(4,1) X X.^2 X.^3];
    c = inv(A)*Y;
end

function c = OrderPolynomialRegression (x0, x1, y0, y1, theta0, theta1)
    A = [1 x0 x0^2 x0^3; 0 1 2*x0 3*x0^2; 1 x1 x1^2 x1^3; 0 1 2*x1 3*x1^2];
    b = [y0; tan(theta0); y1; tan(theta1)];
    %x = inv(A)*b;
    x = A\b;
    c(1) = x(1);
    c(2) = x(2);
    c(3) = x(3);
    c(4) = x(4);
end

function [path, U, dy] = eldm (scenario, indexes, previousPath, P_npDistances, P_ELDM)
    % This function calculates the node points in front of the vehicle,
    % then applies the driver model to calculate the offset to the
    % mid-lane. Finally, a curve is fit onto the node points to result the
    % final path. The entire planning happens in the global UTF frame
    % Inputs:
    % [1] scenario - array with all signals, cut for the planning horizon 
    % [2] previous path - the previously planned trajectory, used for
    % initializing the new trajectory
    % Parameters:    
    %% step 0: calculating the planner frame
    % finding the closest point in the scenario
    distances = ((scenario(:,indexes.X_abs)-previousPath(end,1)).^2+...
        (scenario(:,indexes.Y_abs) - previousPath(end,2)).^2).^0.5;
    nearestPoint = find(distances == min(distances),1);
    
    % estimation of orientation:
    if (size(previousPath,1) > 1)
        theta0 = atan2((previousPath(end,2)-previousPath(end-1,2)),(previousPath(end,1)-previousPath(end-1,1)));
    else
        theta0 = scenario(nearestPoint(1,1), indexes.theta_calc);
    end
    plannerFrame = [previousPath(end,:) theta0];
    Tplanner = [cos(theta0) sin(theta0); -sin(theta0) cos(theta0)];
    
    %% step 1: calculating the node points
    % Node Point Model
    corridorGlobalFrame = [scenario(nearestPoint(1,1):end,indexes.corrX) scenario(nearestPoint(1,1):end,indexes.corrY)];
    corridorPlannerFrame = (corridorGlobalFrame-plannerFrame(1:2))*Tplanner';
    [indeces, ~] = nodePointModel(corridorPlannerFrame, P_npDistances);
    
    % nominal road points
    X_nominal = corridorPlannerFrame(indeces(2:end),1);
    Y_nominal = corridorPlannerFrame(indeces(2:end),2);
    theta_nominal = scenario(indeces(2:end) + nearestPoint(1,1),indexes.LaneOrientation)+ ...
        scenario(indeces(2:end) + nearestPoint(1,1),indexes.theta_calc)- ...
        plannerFrame(3);
    
    %% step 2: calculating the offsets
    % Offset Model
    kappa_nominal = zeros(3,1);
    scenario(:,indexes.LaneCurvature) = movmean(scenario(:,indexes.LaneCurvature),100);
    for j=1:length(indeces)-1
        kappa_nominal(j,1) = mean(scenario(indeces(j) + nearestPoint(1,1):indeces(j+1) + nearestPoint(1,1),indexes.LaneCurvature));
    end
    kappa_nominal = kappa_nominal/0.001;
    % scaling of the predictor variables
    U = [max(kappa_nominal,0); min(kappa_nominal,0); 1];
    % The ouput of the model is the linear combination of predictor variables
    dy = U(1:3)'*P_ELDM(:,1:3)+U(4:6)'*P_ELDM(:,4:6)+P_ELDM(:,7)';
    dy(1:3) = min(max(-1.25,dy(1:3)), 1.25);
    dy = dy';
    
    %% step 3: calculating final node points
    Y = Y_nominal + dy.*cos(theta_nominal);
    X = X_nominal - dy.*sin(theta_nominal);
    theta = theta_nominal;
    
    % extending with origin of the planner frame
    Y = [0; Y]; X = [0;X]; theta = [0; theta];
    
    %% step 4: curve fitting model
    refX = corridorPlannerFrame(:,1)'; %0:1:corridorPlannerFrame(indeces(4),1); %
    [pathPlannerFrame, ~] = trajectory_planner([X;Y;theta], indeces, refX',0);
    
    %% step 5: converting to global frame
    path = pathPlannerFrame*Tplanner + plannerFrame(1:2);
end