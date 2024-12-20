function [path, u, dy] = offsetModel (modelID, scenario, indexes, previousPath, parameters)
    global vehicleState historicalInputData historicalOffset
    

    %% prepare the scenario extended with simulation data for the following fields - and reduce only to the first point
    % vx, ax, yawrate, ay (vehicle kinematics) and LaneOrientation, thetaTP and thetaFP

    syntheticScenario = [scenario(1, indexes.OncomingTrafficType), ...
                            scenario(1, indexes.FrontTrafficType), ...
                            vehicleState.vx, ...
                            vehicleState.ax, ...
                            vehicleState.yawRate, ... 
                            scenario(1, indexes.LaneCurvature), ...
                            scenario(1, indexes.c3), ...
                            vehicleState.thetaTP, ...
                            vehicleState.thetaFP, ...
                            scenario(1, indexes.LaneOrientation), ...
                            vehicleState.ay];

    indexesSynthetic.OncomingTrafficType = 1;
    indexesSynthetic.FrontTrafficType = 2;
    indexesSynthetic.VelocityX = 3;
    indexesSynthetic.AccelerationX = 4;
    indexesSynthetic.YawRate = 5;
    indexesSynthetic.LaneCurvature = 6;
    indexesSynthetic.c3 = 7;
    indexesSynthetic.thetaTP = 8;
    indexesSynthetic.thetaFP = 9;
    indexesSynthetic.LaneOrientation = 10;
    indexesSynthetic.AccelerationY = 11;

    if (indexes.LDM_1 > 0)
        indexesSynthetic.LDM_1 = 12;
        for n=1:indexes.LDM_N-indexes.LDM_1+1
            syntheticScenario = [syntheticScenario scenario(1, indexes.LDM_1+n-1)];            
        end
        indexesSynthetic.LDM_N = indexesSynthetic.LDM_1+indexes.LDM_N-indexes.LDM_1;
    end    
    
    % updating the sweep array of historical knowledge
    if (isempty(historicalInputData))
        for n=1:parameters.PARAMS.HISTORICAL_WINDOW_LENGTH
            historicalInputData(n,:) = [syntheticScenario historicalOffset(end)]; %+1: as the offset is added, as well
        end
    end
    for n=2:size(historicalInputData,1)
        historicalInputData(n,:) = historicalInputData(n-1,:);
    end
    historicalInputData(1,:) = [syntheticScenario historicalOffset(end)];

    %% now generate online input
    inputRaw = prepareDataLive(modelID, syntheticScenario, indexesSynthetic, parameters);

    %% now generate the offset estimate
    [estimation, npDistances, input] = generateEstimate(modelID,inputRaw);

    %% path planning
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
    
    % calculating the node points
    corridorGlobalFrame = [scenario(nearestPoint(1,1):end,indexes.corrX) scenario(nearestPoint(1,1):end,indexes.corrY)];
    corridorPlannerFrame = (corridorGlobalFrame-plannerFrame(1:2))*Tplanner';
    [indeces, ~] = nodePointModel(corridorPlannerFrame, [npDistances(1) diff(npDistances)]/250);
    
    referenceGlobalFrame = [scenario(nearestPoint(1,1):end,indexes.X_abs) scenario(nearestPoint(1,1):end,indexes.Y_abs)];
    referencePlannerFrame = (referenceGlobalFrame-plannerFrame(1:2))*Tplanner';
    Y_reference = referencePlannerFrame(indeces(2:end),2);
    
    % nominal road points
    X_nominal = corridorPlannerFrame(indeces(2:end),1);
    Y_nominal = corridorPlannerFrame(indeces(2:end),2);
    theta_nominal = scenario(indeces(2:end) + nearestPoint(1,1),indexes.LaneOrientation)+ ...
        scenario(indeces(2:end) + nearestPoint(1,1),indexes.theta_calc)- ...
        plannerFrame(3);
    
    dy = estimation'.*cos(theta_nominal);
    u = input';
    
    % calculating final node points
    Y = Y_nominal+ estimation(1,:)'.*cos(theta_nominal);
    X = X_nominal+ estimation(1,:)'.*sin(theta_nominal);
    theta = theta_nominal;
    
    % extending with origin of the planner frame
    Y = [0; Y]; X = [0;X]; theta = [0; theta];
    
    % curve fitting model
    refX = corridorPlannerFrame(1:indeces(end),1)';
    n = 5;
    n = min(n, length(X)-1);
    pathY = fitPolynomialWithConstraints(X,Y, refX, n);
    pathPlannerFrame = [refX' pathY'];
    
    % converting to global frame
    path = pathPlannerFrame*Tplanner + plannerFrame(1:2);
end
