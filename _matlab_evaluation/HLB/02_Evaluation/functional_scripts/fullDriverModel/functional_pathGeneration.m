function [path, U, dY, plannedPath, intentionPath, vehicleStateMemory] = functional_pathGeneration(segment_m, indexes, parameters)
    global  globalStartIndex globalStopIndex vehicleState metadata net dt dts modelMode x_k1 
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
            [priorPath, u, dy] = groundTruth (scenario, indexes, previousPosteriorPath,false, parameters);
        case "groundTruthFiltered"
            [priorPath, u, dy] = groundTruth (scenario, indexes, previousPosteriorPath, true, parameters);
        case "modifiedGroundTruth"
            [priorPath, u, dy] = modifiedGroundTruth (scenario, indexes, previousPosteriorPath);
        case "gp"
            [priorPath, u, dy] = gpModel (scenario, indexes, previousPosteriorPath);
        case "lrm"
            [priorPath, u, dy] = lrmModel (scenario, indexes, previousPosteriorPath);
        otherwise
            priorPath = [scenario(:,indexes.X_abs) scenario(:,indexes.Y_abs)];
    end
end

function [path, u, dy] = groundTruth (scenario, indexes, previousPath, filtered, parameters)
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

function [estimation, deviation] = gpGenerateEstimate(segment_m, indexes)
global parameters
    %% step 0: generate input data
    [~, ~, inputRaw, outputRaw, c_out, s_out, c_in, s_in] = prepareData(segment_m, indexes, parameters.PARAMS);
    
    %% step 1: parameter loading
    % read the learnt data from previous runs
    pathToParams = "C:\git\KDP\publications\GP\results\gp_model";
    paramGP = dir(fullfile(pathToParams,"ETA_input*"));
    for fileID = 1:length(paramGP)
        paramDriverID = str2num(paramGP(fileID).name(strfind(paramGP(fileID).name, 'driver_')+7:strfind(paramGP(fileID).name, '.')-1));
        if (paramDriverID == parameters.PARAMS.DriverID)
            paramData = load(fullfile(paramGP(fileID).folder, paramGP(fileID).name));
            for npID=1:10
                GP_params.hyp_opt_array{npID} = paramData.ETA(npID).hyp_opt;
                GP_params.output_estimation{npID} = paramData.ETA(npID).output_estimation;
                GP_params.input_estimation{npID} = paramData.ETA(npID).input_estimation;
            end
            break;
        else
            paramData = [];
        end
    end

    %% step 3: norm and central
    [inputResim, c, s] = normAndCentral(inputRaw);
    [output, c, s] = normAndCentral(outputRaw);

    %% step 4: generate output
    meanfunc = [];       % Start with a zero mean prior
    eval(strcat('covfunc = ',parameters.PARAMS.KERNEL_TYPE));    % Squared Exponental covariance function
    % ID problem            
    likfunc = @likGauss;    % Gaussian likelihood
    for i = 1:length(GP_params.hyp_opt_array)
        hyp = struct('mean', [], 'cov', 0, 'lik', -1);
        hyp.cov = GP_params.hyp_opt_array{i};
        [estimationGP(:,i), deviationGP(:,i)] = gp(hyp, @infGaussLik, meanfunc, covfunc, likfunc, GP_params.input_estimation{i}, GP_params.output_estimation{i}, inputResim); % extract the mean and covarance functions
    end

    %% step 5: node point calculation
    for i=1:size(estimationGP,2)
        estimation(:,i) = estimationGP(:,i)*s_out(i)+c_out(i);
        deviation(:,i) = deviationGP(:,i)*s_out(i)+c_out(i);
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

function [input, output, inputRaw, outputRaw, c_out, s_out, c_in, s_in] = prepareData(segment_m, indexes, PARAMS)
     % input array
     input = [segment_m(:, indexes.OncomingVehicleTimeToPass), ...
        segment_m(:, indexes.VelocityX), ...
        movmean(segment_m(:, indexes.AccelerationX),20), ...
        movmean(segment_m(:, indexes.YawRate),20), ...
        movmean(segment_m(:, indexes.LaneCurvature), 20), ...
        movmean(segment_m(:, indexes.c3), 200)];
    
    % output array
    output = zeros(size(input,1),numel(PARAMS.OUTPUT_SHIFT));
    for shiftID=1:numel(PARAMS.OUTPUT_SHIFT)
        dT = mean(diff(segment_m(:, indexes.Relative_time)));
        dx = segment_m(:, indexes.VelocityX)*dT;        
        shiftOnOutput = [1:1:size(segment_m(:,indexes.Relative_time),1)]'+floor(PARAMS.OUTPUT_SHIFT(shiftID)./dx);
        for shiftIDonOutput=1:size(input,1)
            if (shiftOnOutput(shiftIDonOutput) > size(input,1))
                break;
            else
                output(shiftIDonOutput,shiftID) = -segment_m(shiftOnOutput(shiftIDonOutput), indexes.c0);
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