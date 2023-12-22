function evaluator_fullDriverModel(segments,config)
  
global segment_m indexes globalStartIndex globalStopIndex vehicleState parameters metadata corFine net refFine modelID nrmsLoss nrmsLossManual dts modelMode

path = [0:1:1000, zeros(1,1001)]';
targetSpeed = 30; %kph
vehiclePath = functional_vehicleModelTest(path, targetSpeed);
    
    %% single simulation
    parameters.P_npDistances = [0.15, 0.2, 0.25];
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

     segment = segments.segments(1).segment;
    [~, segment_m, indexes] = prepareInputForPlanner(segment);
    
    globalStartIndex = 2; % minimum is 2, otherwise it fails!
    globalStopIndex = size(segment_m,1);

    corFine = corridorGeneration(0.5);
    refFine = referenceGeneration(0.5);

    for i=1:length(segments.segments)
    segment = segments.segments(i).segment;
    [~, segment_m, indexes] = prepareInputForPlanner(segment);
    
    globalStartIndex = 2; % minimum is 2, otherwise it fails!
    globalStopIndex = size(segment_m,1);
    
    %% Definition of the parameters
    % Planner: offset model+curve fitting model
    % Curve policy: none

    
    
    % Definition of metadata
    metadata.pathValidity = 1;
    
    corFine = corridorGeneration(0.5);
    refFine = referenceGeneration(0.5);
    
    modelID = "modifiedGroundTruth";
    modelMode = "dynamic";
    %[pathLite, U, dY, ~, intentionPathLite, vehicleStateMemory] = pathGenerationLite();
%     for i=2:length(vehicleStateMemory)
%         steeringAngle(i) = vehicleStateMemory{1,i}.steeringAngle;
%     end
%     disp(strcat('Max difference:',num2str(max(abs(pathLite(globalStartIndex:end,2)-intentionPathLite(globalStartIndex:end,2))))));
%     disp(strcat('Avg difference:',num2str(mean(abs(pathLite(globalStartIndex:end,2)-intentionPathLite(globalStartIndex:end,2))))));
%     disp(strcat('Std difference:',num2str(std(pathLite(globalStartIndex:end,2)-intentionPathLite(globalStartIndex:end,2)))));
    [path, U, dY, ~, intentionPath] = pathGeneration();
    parameters.P_ELDM = reshape(functional_driverModelLearning(U, dY ,8), 3,7);

    parameters_out(i).U = U;
    parameters_out(i).dY = dY;
    parameters_out(i).drID = str2num(segments.segments(i).name(3:5));
    end
    save(fullfile(config.root,"plots", 'parametersAllDrivers.mat'), 'parameters_out');
    
    
    modelID = "eldm";
    parameters.P_curvePolicy(1) = 0.01;
    parameters.P_curvePolicy(2) = 0;
    parameters.P_curvePolicy(3) = 0.75;
    
    path = pathGeneration();

    options = optimset('PlotFcns',@optimplotfval);
    parameters.P_curvePolicy = fminsearch(@optimizationLoss, parameters.P_curvePolicy, options);

    path = pathGeneration();


    %% loop simulation option for different number of node points
%     for i=1:10
%         parameters.P_npDistances = ones(1,i)*(140/i)/250;
%     
%         %parameters.P_npDistances = 135/250; %ones(1,10)*15/250; % [0, 29, 96]/250;
%         parameters.P_LDM = reshape([2.75127773046959;-0.0398873783094187;-1.17733327762608;-2.94882588015453;0.460440066139664;1.33967800733396;0.327136469963189;-0.296796325338061;-0.0234709186554829;0;0;0;0;0;0;0;0;0;0;0;0], 3,7);
%         parameters.P_LDM = parameters.P_LDM(:,1:3);
%         parameters.P_LDM = zeros(length(parameters.P_npDistances));
%         parameters.P_ELDM = reshape([2.86078399132926;0.884814215672794;-1.90657794718284;-3.09943416608130;-0.665457759838954;2.30236448840005;0.348462602099426;-0.107035325513227;-0.271014703397729;1.07959046302992;-0.775251579323662;-0.252977961446196;-0.822164501814478;1.36747233514778;0.113183483561418;-0.124241139196637;-0.454142531428492;0.293625990988783;-0.000983031283019174;-0.000983031283019174;-0.000983031283019174], 3,7);
%         
%         %% Definition of metadata
%         metadata.pathValidity = 1;
%         
%         corFine = corridorGeneration(0.5);
%         refFine = referenceGeneration(0.5);
%          % Init vehicle state
%         vehicleState = initVehicleState(segment_m, indexes, globalStartIndex);
%         path = pathGeneration();
%         [orientation, curvature] = calcPathGeometry(path);
%            
%         offsetPath = offsetCalculation(corFine, path);
%         
%         curves = abs(movmean(curvature,200)) > 2.5e-4;
%         
%         % loss calculation of resulting path
%         indeces = (curves)&(~isnan(offsetPath'));
%         loss = median(abs(offsetPath(indeces==1)))
%         mean(dts)
%         losses(i) = loss;
%         times(i) = mean(dts);
% 
%     end
%     
%     disp(losses);
%     disp(times);

    %% batch parameter simulation
    modelID = "eldm";
    modelMode = "dynamic";
    parameters.P_npDistances = [0.15, 0.2, 0.25];
    paramFiles = dir(fullfile(config.root,'parameters*.mat'));
    averageParamFiles = dir(fullfile(config.root,'averaged*.mat'));
    driverClusters = [2 2 3 2 2 ...
        2 1 1 2 1 ...
        1 1 1 1 1 ...
        1 2 3 2];
    clusterIDs = unique(driverClusters);
    figure(); 
    for k=1:length(clusterIDs)
        subplot(length(clusterIDs),1,k);
        noDriversInCluster = 0;
        P_averaged = [];
        for fileID = 1:length(paramFiles)
            parameterRawData = load(fullfile(paramFiles(fileID).folder,paramFiles(fileID).name));
            if (isfield(parameterRawData, 'parameters'))
                parameterRawData = parameterRawData.parameters;
                fn = fieldnames(parameterRawData);
                numberOfDrivers = length(fn);
            else
                fn = fieldnames(parameterRawData);
                parameterRawData = parameterRawData.(fn{1});
                numberOfDrivers = size(parameterRawData,1);
            end            
            for j=1:numberOfDrivers
                if (length(fn) > 1)
                    drID = fn{j};
                    drID = str2num(drID(3:5));
                    driverParam = parameterRawData.(fn{j});
                else
                    drID = j;
                    driverParam = parameterRawData(j,:);
                    driverParam(end+1:end+2) = driverParam(end);
                end
                if (driverClusters(drID) == clusterIDs(k))
                   noDriversInCluster = noDriversInCluster + 1;
                   parameters.P_ELDM = reshape(driverParam, 3,7);
                   if (isempty(P_averaged))
                       P_averaged = parameters.P_ELDM;
                   else
                       P_averaged = 1/noDriversInCluster * ((noDriversInCluster-1)*P_averaged + parameters.P_ELDM);
                   end
                   [path, ~, ~, ~, ~, vehicleStateMemory] = pathGeneration();
                   modelMode = "kinematic";
                   [pathKinematic, ~, ~, ~, ~, vehicleStateMemoryKinematic] = pathGeneration();
                   offsetPath = offsetCalculation(corFine, path);
                   plot(path(:,1), offsetPath, 'DisplayName', strcat("Dr",num2str(drID))); hold on; grid on;
                   title(strcat('Cluster',{' '}, num2str(clusterIDs(k))));
                end
            end
        end
        parameters.P_ELDM = P_averaged;
        disp(P_averaged);
       [path, U, dY, plannedPath] = pathGeneration();
       offsetPath = offsetCalculation(corFine, path);
       plot(path(:,1), offsetPath, 'DisplayName', 'Averaged', 'LineWidth', 2, 'color', 'k');
       xlabel('UTF-X(m)'); ylabel('offset(m)');
       dataSaveForVisualization(config, U, dY, path, num2str(clusterIDs(k)), segment_m, indexes,offsetPath, plannedPath, P_averaged);

        legend;
        ylim([-1.5, 1.5]);
    end
    
    %% optimization based solution
%     options = optimset('PlotFcns',@optimplotfval);
%     P_optimized = fminsearch(@optimizationLoss, parameters.P_LDM, options);
%     parameters.P_LDM = P_optimized;
%     path = pathGeneration();
%     offsetPath = offsetCalculation(refFine, path);
%     [orientation, curvature] = calcPathGeometry(path);
%     curves = abs(movmean(curvature,200)) > 2.5e-4;
%     
%     % loss calculation of resulting path
%     indeces = (curves)&(~isnan(offsetPath'));
%     loss = median(abs(offsetPath(indeces==1)))

    %% DNN based solution - OPEN-LOOP
    % Simulating whole route to generate target and input vectors
    modelID = "groundTruth";

    parameters.P_npDistances = [0, 29, 96]/250;
    U = []; T=[];
    % Generating train data
    for i=1:size(segments.segments,2)
        segment = segments.segments(i).segment;
        [~, segment_m, indexes] = prepareInputForPlanner(segment);

        globalStartIndex = 2;
        globalStopIndex = size(segment_m,1);
    
        [~, Umem, Tmem] = pathGeneration();
        % modification to get rid of straight line effects
        for j=1:size(Umem,2)
            if (max(abs(Umem(:,j))) <= 0.5)
                straightSection(j) = 1;
            else
                straightSection(j) = 0;
            end
        end
        U = [U Umem(:,straightSection==0)];
        T = [T Tmem(:,straightSection==0)];
        clear straightSection
    end
    Ic = max(abs(U(1:3,:))) <= 0.5;
    U = U(:,Ic==0);
    T = T(:,Ic==0);
    % deleting outliers
    Ic = max(abs(T)) > 1.25;
    U = U(:,Ic==0);
    T = T(:,Ic==0);
    % validation data
    DataV = floor(size(U,2)*0.25);
    Uv = U(:,1:DataV);
    Tv = T(:,1:DataV);
    U = U(:,DataV+1:end);
    T = T(:,DataV+1:end);


    numEpochs = 48000;
    iteration = 0;
    numBatch = 10;
    dataSize = size(U,2);
    dataSizeValidation = size(Uv,2);
    batchSize = floor(dataSize/numBatch);
    batchSizeValidation = floor(dataSizeValidation/numBatch);
    numIterPerBatch = 5;
    iterSize = floor(batchSize/numIterPerBatch);
    iterSizeValidation = floor(batchSizeValidation/numIterPerBatch);
    net = generateNetwork(3);
    averageGrad = []; averageSqGrad = [];
    %learnRate = 0.01;
    %gradDecay = 0.75;

    for epoch = 1:numEpochs
        for batch=1:numBatch
            iteration = iteration + 1;
            for iter=1:numIterPerBatch
                % Prepare mini batch - one scenario simulated only
                % for this, we know the target and the inputs already
                iterStartIdx = (batch-1)*batchSize+(iter-1)*iterSize+1;
                iterStopIdx = iterStartIdx + iterSize;
                dlU(1:3,(iter-1)*iterSize+1:(iter-1)*iterSize+1+iterSize,iter) = U(1:3,iterStartIdx:iterStopIdx);
                dlT(1:3,(iter-1)*iterSize+1:(iter-1)*iterSize+1+iterSize,iter) = T(1:3,iterStartIdx:iterStopIdx);

            end
            [loss,gradients] = dlfeval(@modelLoss,net, dlU, dlT);
            clear dlU dlT
            rmsInEpoch(batch,1) = nrmsLoss;
            
            % validation loss
            
%             e = extractdata(T)-extractdata(V);
%             f = reshape(e,3,size(e,3)*size(e,2));
%             validationLossInEpoch(batch) = norm(f);
            
            %gradients = dlgradient(loss,net.Learnables);
            % Update the network parameters using the Adam optimizer.
            [net,averageGrad,averageSqGrad] = adamupdate(net,gradients, ...
                averageGrad,averageSqGrad,iteration);
        end

        for batch=1:numBatch
            for iter=1:numIterPerBatch
                % Prepare mini batch - one scenario simulated only
                % for this, we know the target and the inputs already
                iterStartIdx = (batch-1)*batchSizeValidation+(iter-1)*iterSizeValidation+1;
                iterStopIdx = iterStartIdx + iterSizeValidation;
                dlUv(1:3,(iter-1)*iterSizeValidation+1:(iter-1)*iterSizeValidation+1+iterSizeValidation, iter) = Uv(1:3,iterStartIdx:iterStopIdx);
                dlTv(1:3,(iter-1)*iterSizeValidation+1:(iter-1)*iterSizeValidation+1+iterSizeValidation,iter) = Tv(1:3,iterStartIdx:iterStopIdx);
            end
            dlU = dlarray(dlUv,"CTB");
            dlT = dlarray(dlTv, "CTB");
            dlY = forward(net,dlU);
            dnnlossValidation = mse(dlY,dlT);
            clear dlUv dlTv
            rmsInEpochValidation(batch,1) = extractdata(dnnlossValidation)^0.5;
        end

        meanRmsForDataSet(epoch) = mean(rmsInEpoch(:,1));   
        meanRmsForDataSetValidation(epoch) = mean(rmsInEpochValidation(:,1));  
      
        subplot(2,1,1);
        plot(meanRmsForDataSet, 'Marker', 'o', 'color', 'k', 'DisplayName', 'train');
        hold on;
        plot(meanRmsForDataSetValidation, 'Marker', 'x', 'color', 'b', 'DisplayName','validation');
        hold off;
        legend;
        
        grid on;
        title(strcat('Current mean loss:', {' '}, num2str(mean(rmsInEpoch))));

        
        if(rem(epoch,100) == 0)
            modelID = "dnn";
            globalStartIndex = 2;
            globalStopIndex = size(segment_m,1);
            vehicleState = initVehicleState(segment_m, indexes, globalStartIndex);
            [~, Um, Tm] = pathGeneration();
            subplot(2,1,2);
            plot(U(1,:), T(1,:),'bo');
            hold on;
            plot(Um(1,:), Tm(1,:),'rx');
            hold off;
            grid on;
        end

        pause(0.1);
    end

    

end

function dlV = generateValidationData(batch, dlV)
global modelID segment_m indexes vehicleState globalStartIndex
    % Init vehicle state
    vehicleState = initVehicleState(segment_m, indexes, globalStartIndex);
    modelID = "dnn";
    [~, ~, V] = pathGeneration();
    if (isempty(dlV))
        clear dlV;
        dlV(1:3,1:size(V,2),batch) = dlarray(V, "CTB");
    else
        dlV(1:3,1:size(V,2),batch) = dlarray(V, "CTB");
    end
end

function dataSaveForVisualization(config, GT_U, GT_Y, path, token, segment_m, indexes, offsetPath, plannedPath, Poptim)
    replan_array = ones(size(path,1),1);
    segment.X_abs = path(:,1);
    segment.Y_abs = path(:,2);
    segment.q_T0 = segment_m(:,indexes.q_T0);
    segment.c2 = segment_m(:,indexes.c2);
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

function [dlU,dlT] = generateDataForNetwork(batch, dlU, dlT)
    global modelID vehicleState segment_m indexes globalStartIndex
    modelID = "groundTruth";
    vehicleState = initVehicleState(segment_m, indexes, globalStartIndex);
    [~, U, T] = pathGeneration();
    if (isempty(dlU))
        clear dlU dlT
        dlU(1:3,1:size(U,2),batch) = dlarray(U, "CTB");
        dlT(1:3,1:size(U,2),batch) = dlarray(T, "CTB");
    else
        dlU(1:3,1:size(U,2),batch) = dlarray(U, "CTB");
        dlT(1:3,1:size(U,2),batch) = dlarray(T, "CTB");
    end
   
    modelID = "dnn";
end

function net = generateNetwork(numChannels)
    % generating using dlnetwork
    layers = [
    sequenceInputLayer(numChannels, 'Name', 'Input')
    fullyConnectedLayer(64, 'Name', 'FC1')
    reluLayer('Name','ReLu1')
    fullyConnectedLayer(32, 'Name', 'FC2')
    reluLayer('Name','ReLu2')
    fullyConnectedLayer(numChannels, 'Name', 'Output')];


    lgraph = layerGraph(layers);

    net = dlnetwork(lgraph);
end

function loss = optimizationLoss(P)
    global parameters segment_m indexes vehicleState globalStartIndex refFine
    parameters.P_curvepolicy = P;
    % Init vehicle state
    vehicleState = initVehicleState(segment_m, indexes, globalStartIndex);
    try
        path = pathGeneration();
    
        offsetPath = offsetCalculation(refFine, path);
        [~, curvature] = calcPathGeometry(path);
        curves = abs(movmean(curvature,200)) > 2.5e-4;

        % loss calculation of resulting path
        indeces = (~isnan(offsetPath'));
        loss = median(abs(offsetPath(indeces==1)));
    catch
        loss = 1000;
    end
end

function [dnnloss, gradients] = modelLoss(net, U, T) 
    global nrmsLoss
    u = dlarray(U,"CTB");
    dlY = forward(net,u);
    t = dlarray(T,"CTB");
    dnnloss = mse(dlY,t);
    nrmsLoss = extractdata(dnnloss)^0.5; %/(mean(mean(mean(extractdata(abs(dlY))))));
    gradients = dlgradient(dnnloss,net.Learnables);
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
    if (numel(find((diff(segment_m(:,indexes.corrX)))>0))==0 && ...
        numel(find((diff(segment_m(:,indexes.corrX)))<0))>0)
        % monotonous decreasing
        x_fine = segment_m(end,indexes.corrX):step:segment_m(1,indexes.corrX);
    elseif (numel(find((diff(segment_m(:,indexes.corrX)))>0))>0 && ...
        numel(find((diff(segment_m(:,indexes.corrX)))<0))==0)
        % monotonous increasing
        x_fine = segment_m(1,indexes.corrX):step:segment_m(end,indexes.corrX);
    else
        %non-monotonous
        x_fine = [];
    end
    
    if (~isempty(x_fine))
        y_fine = spline(segment_m(:,indexes.corrX), segment_m(:,indexes.corrY), x_fine');
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
    replan = 10;
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
        dT = segment_m(i, indexes.q_T0) - segment_m(i-1, indexes.q_T0);
        
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
            window = [50 150];
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
            vehicleState.steeringAngle = controller(posteriorPath, vehicleState);            
            
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
        end
        path(i,1) = vehicleState.X;
        path(i,2) = vehicleState.Y;
        %U = [U u]; dY = [dY dy];
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
    x_k1 = zeros(5,1); x_k1(3) = segment_m(1,indexes.velocityX);

    for i=globalStartIndex:globalStopIndex
        % calculate time step
        dT = segment_m(i, indexes.q_T0) - segment_m(i-1, indexes.q_T0);
        
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
                [c0, c1] = transformVideoData(scenario(1,indexes.c0),scenario(1,indexes.c1), ...
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

% function vehicleState = vehicleModel(vehicleState, dT, mode)
% if(mode == "kinematic")    
%     % kinematic bicycle model
%     % parameters
%     p_wheelBase = 2.7; % meter
%     % calculating yaw-rate based on steering angle
%     vehicleState.yawRate = tan(vehicleState.steeringAngle)/p_wheelBase * vehicleState.vx;
%     vehicleState.theta = vehicleState.theta + dT*vehicleState.yawRate;
%     % displacement
%     dX = vehicleState.vx*dT*cos(vehicleState.theta);
%     dY = vehicleState.vx*dT*sin(vehicleState.theta);
%     % new vehicle position
%     vehicleState.X = vehicleState.X+dX;
%     vehicleState.Y = vehicleState.Y+dY;
%     vehicleState.vy = 0;
%     vehicleState.ax = 0;
%     vehicleState.ay = vehicleState.vx*vehicleState.yawRate;
%     
% elseif (mode == "dynamic")    
%     % dynamic bicycle model
%     % only for high speeds
%     % inputs:
%     % - vehicleState.steeringAngle
%     % - vehicleState.M_av/M_ah - breaking forces
%     % - vehicleState.M_bv/M_bh - breaking forces
%     % parameters:
%     % r: radius of wheels in [m]
%     % c_alfav/h: lateral slip coefficient [N/rad]
%     % c_sv/h: longitudinal slip coefficient [N/%]
%     % c_w: wind coefficient
%     % rho_air: air density
%     % A: preface
%     % J: rotational intertia of the vehicle
%     % m: mass of the vehicle
%     % lv/lh: COG position from front and rear axle
%     % Jwheel: rotiational interatia of the wheels
%     r = 0.309725; c_alfaf = 76500; c_sf = 250; c_alfar = 76500; c_sr = 250;
%     m = 1519;
%     Jwheel = 250; J = 1818;
%     A = 1.5; c_w = 0; rho_air = 1; lf = 1; lr = 1.5; 
%     
% %    vehicleState.s_v = -(vehicleState.v_vx_v - r*vehicleState.drho_v)/(max(abs(vehicleState.v_vx_v), r*vehicleState.drho_v)); % longitudinal slip front, wheel coordinate frame
%     vehicleState.s_f = -(vehicleState.v_fx_v - r*vehicleState.drho_f)/(r*vehicleState.drho_f); % longitudinal slip front, wheel coordinate frame
% 
% %    vehicleState.alfa_v = -vehicleState.v_vy_v /(abs(r*vehicleState.drho_v)); % lateral slip front, wheel coordinate frame
%     vehicleState.alfa_f = -vehicleState.v_fy_v /(r*vehicleState.drho_f); % lateral slip front, wheel coordinate frame
% 
% %    vehicleState.s_h = -(vehicleState.v_hx_v - r*vehicleState.drho_h)/(max(abs(vehicleState.v_hx_v), r*vehicleState.drho_h)); % longitudinal slip front, wheel coordinate frame
%     vehicleState.s_r= -(vehicleState.v_rx_v - r*vehicleState.drho_r)/(r*vehicleState.drho_r); % longitudinal slip front, wheel coordinate frame
% 
% %    vehicleState.alfa_h = -vehicleState.v_hy_v /(abs(r*vehicleState.drho_h)); % lateral slip front, wheel coordinate frame
%     vehicleState.alfa_r = -vehicleState.v_ry_v /(r*vehicleState.drho_r); % lateral slip front, wheel coordinate frame
% 
%     vehicleState.F_fx_v = c_sf*vehicleState.s_f; % longitudinal tyre force front, wheel coordinate frame
%     vehicleState.F_fy_v = c_alfaf*vehicleState.alfa_f; % lateral tyre force front, wheel coordinate frame
%     vehicleState.F_fx = cos(vehicleState.steeringAngle)*vehicleState.F_fx_v - sin(vehicleState.steeringAngle)*vehicleState.F_fy_v; % longitudinal tyre force front, vehicle frame
%     vehicleState.F_fy = sin(vehicleState.steeringAngle)*vehicleState.F_fx_v + cos(vehicleState.steeringAngle)*vehicleState.F_fy_v; % lateral tyre force front, vehicle frame
%     vehicleState.F_rx = c_sr*vehicleState.s_r; % longitudinal tyre force, rear wheel, vehicle coordinate frame == wheel coordinate frame
%     vehicleState.F_ry = c_alfar*vehicleState.alfa_r; % lateral tyre force, rear wheel, vehicle coordinate frame == wheel coordinate frame
%     vehicleState.F_wx = 0.5*c_w*rho_air*A*vehicleState.v_fx^2;
%     vehicleState.F_wy = 0.5*c_w*rho_air*A*vehicleState.v_fy^2;
%     
%     % COG quantities    
%     vehicleState.a_x = 1/m*(vehicleState.F_fx + vehicleState.F_rx - vehicleState.F_wx);
%     vehicleState.a_y = 1/m*(vehicleState.F_fy + vehicleState.F_ry - vehicleState.F_wy);
%     
%     vehicleState.v_x = vehicleState.v_x + vehicleState.a_x*dT;
%     vehicleState.v_y = vehicleState.v_y + vehicleState.a_y*dT;
%     
%     vehicleState.eps_sigma = 1/J*(lf*vehicleState.F_fy - lr*vehicleState.F_ry);
%     vehicleState.yawRate = vehicleState.yawRate + dT*vehicleState.eps_sigma;
%     
%     % Baselink quantities
%     vehicleState.ax = vehicleState.a_x;
%     vehicleState.vx = vehicleState.v_x;
%     vehicleState.ay = vehicleState.yawRate*vehicleState.vx;
%     vehicleState.vy = 0;
%     
%     % transforming to front and rear wheel
%     vehicleState.v_fx = vehicleState.v_x;
%     vehicleState.v_fy = vehicleState.v_y + vehicleState.yawRate*lf;
%     vehicleState.v_rx = vehicleState.v_x;
%     vehicleState.v_ry = vehicleState.v_y - vehicleState.yawRate*lr;
%     
%     vehicleState.v_fx_v = cos(vehicleState.steeringAngle)*vehicleState.v_fx + sin(vehicleState.steeringAngle)*vehicleState.v_fy;
%     vehicleState.v_fy_v = -sin(vehicleState.steeringAngle)*vehicleState.v_fx + cos(vehicleState.steeringAngle)*vehicleState.v_fy;
%     vehicleState.v_rx_v = vehicleState.v_rx;
%     vehicleState.v_ry_v = vehicleState.v_ry;
%     
%     vehicleState.ddrho_f = 1/Jwheel * (vehicleState.M_af - vehicleState.M_bf)*sign(vehicleState.drho_f) - r*vehicleState.F_fx_v;
%     vehicleState.ddrho_r = 1/Jwheel * (vehicleState.M_ar - vehicleState.M_br)*sign(vehicleState.drho_r) - r*vehicleState.F_rx;
%     
%     vehicleState.drho_f = vehicleState.drho_f + vehicleState.ddrho_f*dT;
%     vehicleState.drho_r = vehicleState.drho_r + vehicleState.ddrho_r*dT;
%     
%     % absolute frame quantities
%     vehicleState.theta = vehicleState.theta + dT*vehicleState.yawRate;
%     % displacement
%     dX = vehicleState.v_rx*dT*cos(vehicleState.theta)-sin(vehicleState.theta)*vehicleState.v_ry*dT;
%     dY = vehicleState.v_rx*dT*sin(vehicleState.theta)+cos(vehicleState.theta)*vehicleState.v_ry*dT;
%     % new vehicle position
%     vehicleState.X = vehicleState.X+dX;
%     vehicleState.Y = vehicleState.Y+dY;
%     
%     vehicleState.v = (vehicleState.v_x^2 + vehicleState.v_y^2)^0.5;  
%     
% elseif (mode == "dynamicSimplified")
%     r = 0.309725; c_alfaf = 76500; c_sf = 250; c_alfar = 76500; c_sr = 250;
%     m = 1519;
%     Jwheel = 250; J = 1818;
%     lf = 1; lr = 1.5;     
%     
%     df = vehicleState.steeringAngle;
%     vehicleState.F_fx = cos(df)*(-c_sf*(cos(df)*vehicleState.v_x+sin(df)*vehicleState.v_y+sin(df)*vehicleState.yawRate*lf-r*vehicleState.drho_f)/(r*vehicleState.drho_f))+ ...
%         sin(df)*c_alfaf*(-sin(df)*vehicleState.v_x+cos(df)*vehicleState.v_y+cos(df)*vehicleState.yawRate*lf)/(r*vehicleState.drho_f);
%     vehicleState.F_fy = -sin(df)*c_sf*(cos(df)*vehicleState.v_x+sin(df)*vehicleState.v_y+sin(df)*vehicleState.yawRate*lf-r*vehicleState.drho_f)/(r*vehicleState.drho_f) - ...
%         cos(df)*c_alfaf*(-sin(df)*vehicleState.v_x+cos(df)*vehicleState.v_y+cos(df)*vehicleState.yawRate*lf)/(r*vehicleState.drho_f);
%     vehicleState.F_rx = c_sr*(-(vehicleState.v_x-r*vehicleState.drho_r)/(r*vehicleState.drho_f));
%     vehicleState.F_ry = c_alfar*(-(vehicleState.v_y-vehicleState.yawRate*lr)/(r*vehicleState.drho_r));
%     vehicleState.a_x = 1/m*(vehicleState.F_fx + vehicleState.F_rx);
%     vehicleState.a_y = 1/m*(vehicleState.F_fy + vehicleState.F_ry);
%     vehicleState.v_x = vehicleState.v_x + vehicleState.a_x*dT;
%     vehicleState.v_y = vehicleState.v_y + vehicleState.a_y*dT;
%     vehicleState.yawRate = vehicleState.yawRate + dT/J*lf*vehicleState.F_fy - dT*lr/J*vehicleState.F_ry;
%     vehicleState.drho_f = vehicleState.drho_f+dT*r*c_sf*(cos(df)*vehicleState.v_x+sin(df)*vehicleState.v_y+sin(df)*vehicleState.yawRate*lf-r*vehicleState.drho_f)/(r*vehicleState.drho_f);
%     vehicleState.drho_r = vehicleState.drho_r+r*dT*c_sr*(vehicleState.vx-r*vehicleState.drho_r)/(r*vehicleState.drho_r);
%     
%     vehicleState.theta = vehicleState.theta + vehicleState.yawRate*dT;
%     % displacement
%     v_hx = vehicleState.v_x;
%     v_hy = vehicleState.v_y - lr*vehicleState.yawRate;
%     dX = v_hx*dT*cos(vehicleState.theta)-sin(vehicleState.theta)*v_hy*dT;
%     dY = v_hx*dT*sin(vehicleState.theta)+cos(vehicleState.theta)*v_hy*dT;
%     % new vehicle position
%     vehicleState.X = vehicleState.X+dX;
%     vehicleState.Y = vehicleState.Y+dY;
%     
% end
% end

function vehicleStateCheck = checkModelEquality(vehicleState,vehicleStateCheck)
    vehicleStateCheck.F_fx_err = vehicleStateCheck.F_fx - vehicleState.F_fx;
    vehicleStateCheck.F_fy_err = vehicleStateCheck.F_fy - vehicleState.F_fy;
    vehicleStateCheck.F_rx_err = vehicleStateCheck.F_rx - vehicleState.F_rx;
    vehicleStateCheck.F_ry_err = vehicleStateCheck.F_ry - vehicleState.F_ry;
    vehicleStateCheck.drho_f_err = vehicleStateCheck.drho_f - vehicleState.drho_f;
    vehicleStateCheck.drho_r_err = vehicleStateCheck.drho_r - vehicleState.drho_r;
end

function vehicleState = initVehicleState(segment_m, indexes, index, mode)
if (mode=="kinematic")
    vehicleState.X = segment_m(index, indexes.X_abs);
    vehicleState.Y = segment_m(index, indexes.Y_abs);
    vehicleState.theta = segment_m(index, indexes.theta_calc);
    vehicleState.yawRate = segment_m(index, indexes.yawRate);
    vehicleState.vx = segment_m(index, indexes.velocityX);
    vehicleState.vy = 0;
    vehicleState.ax = 0;
    vehicleState.ay = vehicleState.vx*vehicleState.yawRate;
    vehicleState.t0 = segment_m(index, indexes.q_T0);    
    vehicleState.steeringAngle = 0;
elseif (mode =="dynamic")
    lf = 1; lr = 1.5; r = 0.309725;
    
    vehicleState.X = segment_m(index, indexes.X_abs);
    vehicleState.Y = segment_m(index, indexes.Y_abs);
    vehicleState.theta = segment_m(index, indexes.theta_calc);
    vehicleState.v_x = segment_m(index, indexes.velocityX);
    vehicleState.v_y = 0;
    vehicleState.t0 = segment_m(index, indexes.q_T0);
    vehicleState.yawRate = segment_m(index, indexes.yawRate);
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
elseif (mode =="dynamicSimplified")
    lf = 1; lr = 1.5; r = 0.309725;
    
    vehicleState.X = segment_m(index, indexes.X_abs);
    vehicleState.Y = segment_m(index, indexes.Y_abs);
    vehicleState.theta = segment_m(index, indexes.theta_calc);
    vehicleState.v_x = segment_m(index, indexes.velocityX);
    vehicleState.v_y = 0;
    vehicleState.t0 = segment_m(index, indexes.q_T0);
    vehicleState.yawRate = segment_m(index, indexes.yawRate);
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
        case "dnn"
            [priorPath, u, dy] = dnn (scenario, indexes, previousPosteriorPath, net);
        case "groundTruth"
            [priorPath, u, dy] = groundTruth (scenario, indexes, previousPosteriorPath);
        case "modifiedGroundTruth"
            [priorPath, u, dy] = modifiedGroundTruth (scenario, indexes, previousPosteriorPath);
        otherwise
            priorPath = [scenario(:,indexes.X_abs) scenario(:,indexes.Y_abs)];
    end
end

function posteriorPath = curvePolicy (scenario, indexes, priorPath, vehicleState)
    global modelID parameters
    % curve policy returns the posterior path which is a corrigated path
    % based on the control dynamics of a given driver
    % transforming prior Path to local frame
    
    if (modelID ~="groundTruth")
        T = [cos(vehicleState.theta) sin(vehicleState.theta); -sin(vehicleState.theta) cos(vehicleState.theta)];
        priorPathEgoFrame = (priorPath-[vehicleState.X vehicleState.Y])*T';
        posteriorPath = priorPathEgoFrame;
        P.Kp = parameters.P_curvePolicy(1);
        P.Ki = parameters.P_curvePolicy(2);
        % time channel production, considering current vehicle speed
        vx = vehicleState.vx;
        displacements = (diff(priorPath(:,1)).^2 + diff(priorPath(:,2)).^2).^0.5; % in meters
        dtimes = displacements/vx;
        % find the look ahead point
        lad = vx*parameters.P_curvePolicy(3); % the look ahead time
        % loop through the prioPath to modify it according to control driver
        nearestIndex = find(posteriorPath(:,1)>=0,1);
        
        modelVehicleState = vehicleState;
        modelVehicleState.X = 0;
        modelVehicleState.Y = 0;
        modelVehicleState.theta = 0;
        ie = 0;
        for i=max(2,nearestIndex(1,1)):size(priorPath,1)
            yR = spline(priorPathEgoFrame(i-1:end,1), priorPathEgoFrame(i-1:end,2), lad+posteriorPath(i-1,1));
            [modelVehicleState.steeringAngle, ie] = compensatoryDriverModel(yR, modelVehicleState.Y, P, dtimes(i-1), ie);
            modelVehicleState = vehicleModel(modelVehicleState, dtimes(i-1));
            posteriorPath(i,1) = modelVehicleState.X;
            posteriorPath(i,2) = modelVehicleState.Y;
        end
        % transforming back to global frame
        posteriorPath = posteriorPath*T + [vehicleState.X vehicleState.Y];
    else
        posteriorPath = priorPath;
    end

end

function [steeringAngle, i] = compensatoryDriverModel(yref, y, P_controlDriverModel, Ts, ik1)
    e = yref-y;
    i = ik1 + Ts*e;
    steeringAngle = P_controlDriverModel.Kp*e + P_controlDriverModel.Ki*i;
end

function [path, u, dy] = dnn (scenario, indexes, previousPath, net)
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
    
    referenceGlobalFrame = [scenario(nearestPoint(1,1):end,indexes.X_abs) scenario(nearestPoint(1,1):end,indexes.Y_abs)];
    referencePlannerFrame = (referenceGlobalFrame-plannerFrame(1:2))*Tplanner';
    Y_reference = referencePlannerFrame(indeces(2:end),2);
    
    % nominal road points
    X_nominal = corridorPlannerFrame(indeces(2:end),1);
    Y_nominal = corridorPlannerFrame(indeces(2:end),2);
    theta_nominal = scenario(indeces(2:end) + nearestPoint(1,1),indexes.c1)+ ...
        scenario(indeces(2:end) + nearestPoint(1,1),indexes.theta_calc)- ...
        plannerFrame(3);
    
    %% step 2: calculating the offsets
    % Offset Model
    u = zeros(3,1);
    for j=1:length(indeces)-1
        u(j,1) = mean(scenario(indeces(j) + nearestPoint(1,1):indeces(j+1) + nearestPoint(1,1),indexes.c2));
    end
    % scaling of the predictor variables
    kappa_lim = 0.001;
    U_lim = ones(3,1)*kappa_lim;
    u = u./U_lim;

    % The ouput of the model is the linear combination of predictor variables
    u_(1:3,1,1) = u;
    dlu = dlarray(u_, 'CTB');
    dy = forward(net,dlu);
    
    dy = extractdata(dy);
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

function [path, u, dy] = groundTruth (scenario, indexes, previousPath)
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
    
    referenceGlobalFrame = [scenario(nearestPoint(1,1):end,indexes.X_abs) scenario(nearestPoint(1,1):end,indexes.Y_abs)];
    referencePlannerFrame = (referenceGlobalFrame-plannerFrame(1:2))*Tplanner';
    Y_reference = referencePlannerFrame(indeces(2:end),2);
    
    % nominal road points
    X_nominal = corridorPlannerFrame(indeces(2:end),1);
    Y_nominal = corridorPlannerFrame(indeces(2:end),2);
    theta_nominal = scenario(indeces(2:end) + nearestPoint(1,1),indexes.c1)+ ...
        scenario(indeces(2:end) + nearestPoint(1,1),indexes.theta_calc)- ...
        plannerFrame(3);
    
    %% step 2: calculating the offsets
    % Offset Model
    u = zeros(3,1);
    kappa_nominal = zeros(3,1);
    scenario(:,indexes.c2) = movmean(scenario(:,indexes.c2),100);
    for j=1:length(indeces)-1
        kappa_nominal(j,1) = mean(scenario(indeces(j) + nearestPoint(1,1):indeces(j+1) + nearestPoint(1,1),indexes.c2));
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
    theta_nominal = scenario(indeces(2:end) + nearestPoint(1,1),indexes.c1)+ ...
        scenario(indeces(2:end) + nearestPoint(1,1),indexes.theta_calc)- ...
        plannerFrame(3);
    
    %% step 2: calculating the offsets
    % Offset Model
    u = zeros(3,1);
    kappa_nominal = zeros(3,1);
    scenario(:,indexes.c2) = movmean(scenario(:,indexes.c2),100);
    for j=1:length(indeces)-1
        kappa_nominal(j,1) = mean(scenario(indeces(j) + nearestPoint(1,1):indeces(j+1) + nearestPoint(1,1),indexes.c2));
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
    
    scenario(:,indexes.c2) = movmean(scenario(:,indexes.c2),100);
    
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
    theta_nominal = scenario(indeces(2:end) + nearestPoint(1,1),indexes.c1)+ ...
        scenario(indeces(2:end) + nearestPoint(1,1),indexes.theta_calc)- ...
        plannerFrame(3);
    
    %% step 2: calculating the offsets
    % Offset Model
    U = zeros(N,1);
    for j=1:length(indeces)-1
        U(j,1) = mean(scenario(indeces(j) + nearestPoint(1,1):indeces(j+1) + nearestPoint(1,1),indexes.c2));
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
    N = 3; % number of node points
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
%             u(i,1) = 2*scenario(1,indexes.c2)+3*scenario(1,indexes.c3)*P_npDistances(i);
%         else
%             u(i,1) = 2*scenario(1,indexes.c2)+3*scenario(1,indexes.c3)*(P_npDistances(i)+P_npDistances(i-1));
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
        u(j,1) = mean(scenario(indeces(j):indeces(j+1),indexes.c2));
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
        theta_nominal = atan(scenario(indeces(2:end),indexes.c1))+scenario(indeces(2:end),indexes.theta_calc)-vehicleState.theta;
    else
        for i=1:N
            X_nominal(i) = P_npDistances(i);
            Y_nominal(i) = c0 + c1*P_npDistances(i) + ...
                scenario(1,indexes.c2)/2*P_npDistances(i)^2+ ...
                scenario(1,indexes.c3)/6*P_npDistances(i)^3;
            theta_nominal(i) = atan(c1 + ...
                scenario(1,indexes.c2)*P_npDistances(i)+ ...
                scenario(1,indexes.c3)/2*P_npDistances(i)^2);
        end
    end
            
    
    %% step 2: offset model
    % scaling of the predictor variables
    kappa_lim = 0.001;
    U_lim = [ones(6,1)*kappa_lim; 1];
    U = U./U_lim;
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
    theta_nominal = scenario(indeces(2:end) + nearestPoint(1,1),indexes.c1)+ ...
        scenario(indeces(2:end) + nearestPoint(1,1),indexes.theta_calc)- ...
        plannerFrame(3);
    
    %% step 2: calculating the offsets
    % Offset Model
    U = zeros(7,1);
    kappa_nominal = zeros(3,1);
    scenario(:,indexes.c2) = movmean(scenario(:,indexes.c2),100);
    for j=1:length(indeces)-1
        kappa_nominal(j,1) = mean(scenario(indeces(j) + nearestPoint(1,1):indeces(j+1) + nearestPoint(1,1),indexes.c2));
    end
    % scaling of the predictor variables
    kappa_lim = 0.001;
    U_lim = [ones(6,1)*kappa_lim; 1];
    U = [max(kappa_nominal,0); min(kappa_nominal,0); 1];
    U = U./U_lim;
    % The ouput of the model is the linear combination of predictor variables
    dy = P_ELDM*U;
    dy(1:3) = min(max(-1.25,dy(1:3)), 1.25);
    
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