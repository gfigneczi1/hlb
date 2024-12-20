function evaluator_fullDriverModel_GPProblem(segments,config)
SIMULATION_MODE = "SYNTHETIC_SIMULATION";
global segment_m indexes globalStartIndex globalStopIndex vehicleState parameters metadata corFine net refFine modelID nrmsLoss nrmsLossManual dts modelMode params eliminatedOffsets U_train Y_train

addpath(fullfile("..", "..", "Solutions", "CruiseConcept", "gpBasedPathPlanner", "library"));
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
parameters.kappa_min = 2.5e-4; % in 1/m unit
parameters.eliminateOffset = true;

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
PARAMS.npDistances = [10, 39, 136];
PARAMS.DriverID = 4;
PARAMS.FILTER_OUTPUT = false;
PARAMS.Lp = 1.5;
PARAMS.HISTORICAL_WINDOW_LENGTH = 15;

modelID_array = ["ldm", "eldm", "lrm", "sgp", "phtpm"];
parameters.PARAMS = PARAMS;

colors = ["r", "b", "c", "g", "k", "m", "r--", "b--", "c--", "g--", "k--", "m--", "r:","b:"];
load('C:\database\KDP_HLB_GP\Dr008_Dr027_input_GP.mat');
for segmentID = 5:length(segments.segments)
    if (segmentID ~=2)

        segment = segments.segments(segmentID).segment;
        name = segments.segments(segmentID).name;
        driverID = segmentID;

        parameters.PARAMS.DriverID = driverID;
        PARAMS.DriverID = driverID;
        [~, segment_m, indexes] = prepareInputForPlanner(segment);

        modelMode = "kinematic";
                
        xRel = [1418000, 1421000];

        if (mean(diff(segment_m(:,indexes.X_abs))) < 0)
            globalStartIndex = find(segment_m(:,indexes.X_abs) < xRel(2),1);
            globalStopIndex = find(segment_m(:,indexes.X_abs) < xRel(1),1);
        else
            globalStartIndex = find(segment_m(:,indexes.X_abs) > xRel(1),1);
            globalStopIndex = find(segment_m(:,indexes.X_abs) > xRel(2),1);
        end

        for mId=1:length(modelID_array)
            modelID = modelID_array(mId);
            params = loadParameters(modelID, driverID);
            [segment_m, indexes, expectedOutput] = conditionData(modelID, segment_m, indexes, PARAMS);
        
            [path, U, dY, ~, ~, ~, corridorResampled] = pathGeneration();
    
            if (mean(diff(segment_m(:,indexes.X_abs))) < 0)
                % 62A - east to west
                dir = 1;
            else
                dir = 2;
            end
    
            offsets{driverID, mId}.offset = -1*(path(find(metadata.pathValidity==1)-globalStartIndex+1,2)-corridorResampled(find(metadata.pathValidity==1)-globalStartIndex+1,2)); 
            offsets{driverID, mId}.name = strcat('Dr', num2str(driverID)); 
            offsets{driverID, mId}.X = path(find(metadata.pathValidity==1)-globalStartIndex+1,1); 
            offsets{driverID, mId}.marker = colors(driverID);
            offsets{driverID, mId}.dY = dY(:,find(metadata.pathValidity==1)-globalStartIndex+1,1);
            offsets{driverID, mId}.U = U(:, find(metadata.pathValidity==1)-globalStartIndex+1,1);
            offsets{driverID, mId}.params = params;
            if (dir == 1)
                offsets{driverID, mId}.referenceX = segment_m(:,indexes.X_abs);
                offsets{driverID, mId}.referenceOffset = -segment_m(:,indexes.c0)-params.c_out(1);
            else
                offsets{driverID, mId}.referenceX = segment_m(:,indexes.X_abs);
                offsets{driverID, mId}.referenceOffset = segment_m(:,indexes.c0)+params.c_out(1);
            end
        end

        % now making the plot
        f = figure(1);
        f.Position = [100 100 1050 750];
        set(f,'defaulttextInterpreter','latex') ;
        set(f, 'defaultAxesTickLabelInterpreter','latex');  
        set(f, 'defaultLegendInterpreter','latex');
        for mId = 1:size(offsets,2)
            subplot(size(offsets,2),1,mId);
            plot(offsets{driverID, mId}.referenceX, offsets{driverID, mId}.referenceOffset, 'color', 'k', 'DisplayName','Reference');
            hold on;
            plot(offsets{driverID, mId}.X, offsets{driverID, mId}.offset, 'color', 'b', 'DisplayName', 'Prediction');
            xlim(xRel);
            ylim([-0.6,0.6]);
            title(strcat("Model:", {' '}, modelID_array(mId), ",", {' '}, offsets{driverID, mId}.name));
            grid on;
            legend;
        end

        clear offsets
    end
end

f = plot_offsetPlotsGP_model(offsets, path);
temp_folder_path = "../../_temp";
plots_folder_name = "plots";
savefig(f, fullfile(temp_folder_path, plots_folder_name,...
                    'Resimulate_plots_syntheticData.fig'));
saveas(f, fullfile(temp_folder_path, plots_folder_name,...
                'Resimulate_plots_syntheticData.png'));
close(f);

end

%% *************** END OF MAIN FUNCTION **********************

function [path, U, dY, plannedPath, intentionPath, vehicleStateMemory, corridorResampled] = pathGeneration()
    global segment_m indexes globalStartIndex globalStopIndex vehicleState parameters metadata net dt dts modelMode x_k1 modelID historicalInputData historicalOffset U_train Y_train
    % Entire model is calculated in the global frame!
    vehicleState = initVehicleState(segment_m, indexes, globalStartIndex, modelMode);
    replan = true;
    scenario = []; % initializing the scenario
    metadata= [];
    historicalInputData = [];
    historicalOffset = 0;
    path(1,1) = vehicleState.X;
    path(1,2) = vehicleState.Y;
    filteredOffset = segment_m(globalStartIndex, indexes.c0);
    corridorResampled = [segment_m(globalStartIndex,indexes.corrX) segment_m(globalStartIndex,indexes.corrY)];
    replanCounter = 0;
    U = []; dY = []; plannedPath = []; dts = []; intentionPath = [];
    priorPath = zeros(1,2);
    globalNewCorridorPoint = [0 0];
    u = zeros(3,1); dy = zeros(3,1); 
    x_k1 = zeros(4,1);
    i = globalStartIndex;
    r = 1;
    finished = false;
    routeLength = sum(sqrt(diff(segment_m(globalStartIndex:globalStopIndex,indexes.X_abs)).^2+diff(segment_m(globalStartIndex:globalStopIndex,indexes.Y_abs)).^2));
    while (~finished)
    %for i=globalStartIndex:globalStopIndex
       
        % calculate time step
        dT = 0.1; %segment_m(i, indexes.Relative_time) - segment_m(i-1, indexes.Relative_time);
        
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
                    U_train{r,1} = u;
                    Y_train{r,1} = dy;
                    r = r + 1;
                catch
                    % input data is corrupt and planner failed
                    scenario = [];
                    replan = 1;
                end
                
                dts = [dts; dt];
            end
        end
        % from this point on, only calculate if scenario is valid!!
        if (~isempty(scenario) && dT <= 0.1)
            % curve policy: a model which takes the input data and outputs the
            % modified points
            % - cuts the same subsegment as done for the curve
            % get the nearest index again for the moved vehicle
            %nearestIndex = getNearestIndex(segment_m, indexes, vehicleState);
            % do the cut...
            scenarioCurvePolicy = scenario; %cutScenario(segment_m, indexes, nearestIndex, vehicleState, window);
            % calculate the posterior path based on the curve policy
            posteriorPath = priorPath;
            
            % calculate vehicle position - supposing perfect controller
            % calculate lane orientation
            theta = atan(diff(posteriorPath(1:2,2))/diff(posteriorPath(1:2,1)));
            if (diff(posteriorPath(1:2,1)) < 0 && diff(posteriorPath(1:2,2)) < 0)
                % third quarter plane
                theta = theta - pi();
            elseif (diff(posteriorPath(1:2,1)) < 0 && diff(posteriorPath(1:2,2)) > 0)
                % second quarter plane
                theta = theta + pi();
            end
            % converting the posteriorPath to local coordinates
            T = [cos(theta) sin(theta); -sin(theta) cos(theta)];
            posteriorLocalPath = (posteriorPath - [vehicleState.X vehicleState.Y])*T';
            corridorLocally = ([segment_m(:, indexes.corrX) segment_m(:,indexes.corrY)] - [vehicleState.X vehicleState.Y])*T';

            dispDuringOneCycle = vehicleState.vx*dT;
            Y = spline(posteriorLocalPath(:,1), posteriorLocalPath(:,2), dispDuringOneCycle);
            Ycorridor = spline(corridorLocally(:,1), corridorLocally(:,2), dispDuringOneCycle);

            % interpolating new position
            globalNewPoint = [dispDuringOneCycle Y]*T+[vehicleState.X vehicleState.Y];
            globalNewCorridorPoint = [dispDuringOneCycle Ycorridor]*T+[vehicleState.X vehicleState.Y];
            
            vehicleState.X = globalNewPoint(1);
            vehicleState.Y = globalNewPoint(2);
            vehicleState.theta = theta;

                      
            % finding the closest point
            distances = sqrt((scenario(:,indexes.corrX)-globalNewCorridorPoint(1)).^2+(scenario(:,indexes.corrY)-globalNewCorridorPoint(2)).^2);
            roadCurvature = scenario(find(distances==min(distances),1), indexes.LaneCurvature);
            vehicleState.vx = scenario(find(distances==min(distances),1), indexes.VelocityX);
            vehicleState.ax = scenario(find(distances==min(distances),1), indexes.AccelerationX);
            
            if (isfield(indexes, "thetaTP"))
                vehicleState.thetaTP = scenario(find(distances==min(distances),1), indexes.thetaTP);
                vehicleState.thetaFP = scenario(find(distances==min(distances),1), indexes.thetaFP);
            else
                vehicleState.thetaTP = 0;
                vehicleState.thetaFP = 0;
            end

            % calculating yaw-rate based on steering angle
            vehicleState.yawRate = scenario(find(distances==min(distances),1), indexes.YawRate);
            vehicleState.vy = 0;            
            vehicleState.ay = vehicleState.vx*vehicleState.yawRate;

            filteredOffset = 0.9*filteredOffset+0.1*Ycorridor;

            historicalOffset = [historicalOffset; -scenario(find(distances==min(distances),1), indexes.c0)];
            %historicalOffset = [historicalOffset; -filteredOffset];% multiple by -1 to have the same lane offset sign conventions
            
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
                finished = true;
                break;
            end
            replan = 1;
            metadata.pathValidity(i) = 0;
            replanCounter = 0;
            priorPath = path(end,:);
            posteriorPath = priorPath;
            if (~isempty(historicalOffset))
                historicalOffset = [historicalOffset; historicalOffset(end)];
            end
        end
        pathLength = sum(sqrt(diff(path(:,1)).^2+diff(path(:,2)).^2));
        path = [path; [vehicleState.X vehicleState.Y]];
        U = [U u]; dY = [dY dy];
        intentionPath = [intentionPath; priorPath(1,:)];
        plannedPath = [plannedPath; posteriorPath(1,:)];
        corridorResampled = [corridorResampled; globalNewCorridorPoint];
        vehicleStateMemory{i} = vehicleState;
        if (pathLength >= routeLength || i >= globalStopIndex)
            finished = true;
        end
        i = i+1;
    end
end

function vehicleState = initVehicleState(segment_m, indexes, index, mode)
    vehicleState.X = segment_m(index, indexes.X_abs);
    vehicleState.Y = segment_m(index, indexes.Y_abs);
    vehicleState.theta = segment_m(index, indexes.theta_calc);
    vehicleState.yawRate = segment_m(index, indexes.YawRate);
    vehicleState.vx = segment_m(index, indexes.VelocityX);
    if (isfield(indexes, "thetaTP"))
        vehicleState.thetaTP = segment_m(index, indexes.thetaTP);
        vehicleState.thetaFP = segment_m(index, indexes.thetaFP);
    else
        vehicleState.thetaTP = 0;
        vehicleState.thetaFP = 0;
    end
    vehicleState.vy = 0;
    vehicleState.ax = 0;
    vehicleState.ay = vehicleState.vx*vehicleState.yawRate;

    vehicleState.t0 = segment_m(index, indexes.Relative_time);    
    vehicleState.steeringAngle = 0;
    vehicleState.time = 0;
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
    u = []; dy = [];
    [priorPath, u, dy] = offsetModel (modelID, scenario, indexes, previousPosteriorPath, parameters);
end

function [dataOut, c,s]= normAndCentral(dataIn)
    for i=1:size(dataIn,2)
        c(i) = mean(dataIn(:,i));
        s(i) = std(dataIn(:,i));
        dataOut(:,i) = (dataIn(:,i)-mean(dataIn(:,i)))/std(dataIn(:,i));
    end
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

function theta = calculateOrientation(x,y)
    theta = atan(diff(y)./diff(x)); theta = [theta; theta(end)];
    % check quarterplane 2 and 3 conditions
    dx = diff(x); dx = [dx; dx(end)];
    dy = diff(y); dy = [dy; dy(end)];
    
    qp2 = (dx < 0) & (dy > 0); % the division in the atan function is negative, therefore the angle is negative - compensation of +pi() is needed
    qp3 = (dx < 0) & (dy < 0); % the division in the atan function is positive, therefore the angle is positive - compensation of -pi() is needed
    
    theta(qp2) = theta(qp2)+pi();
    theta(qp3) = theta(qp3)-pi();
end