function [priorPath, u, dy] = functional_planner (scenario, indexes, previousPosteriorPath, parameters)
    global modelID
    % some planner, which returns the priorPath
    % priorPath = NX2 array, path points ahead of the vehicle with a certain
    % step-size (1m)
    % priorPath = [scenario(:, indexes.X_abs) scenario(:, indexes.Y_abs)]; % temporary: ground-truth
    P_npDistances = parameters.P_npDistances;
    u = []; dy = [];
    switch modelID
        case "eldm"
            [priorPath, u ,dy] = eldm (scenario, indexes, previousPosteriorPath, P_npDistances, parameters.P_ELDM);
        case "groundTruth"
            [priorPath, u, dy] = groundTruth (scenario, indexes, previousPosteriorPath,false);
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