function evaluator_optimalPath(segments,config)

% This evaluator is responsible for generating various plots to analyse
% lane wandering performance of different drivers.
% The following KPI-s are proposed to be used:
% - straight line:
% -- compensation points (where the steering torque is applied to compensate stray effect)
% -- distribution of lane offset and CLP
% -- oncoming traffic? follow traffic?

% Segmentor: segmentor_driverModel.

% @copyright (C) 2022 Robert Bosch GmbH.
% The reproduction, distribution and utilization of this file as
% well as the communication of its contents to others without express
% authorization is prohibited. Offenders will be held liable for the
% payment of damages. All rights reserved in the event of the grant
% of a patent, utility model or design.
% @version 1.0
% @date 2023-09-04


pathPlanning (segments, config);

end

function pathPlanning (segments, config)
global modelMode globalStartIndex globalStopIndex indexes segment_m parameters config

modelMode = "kinematic";

segment = segments.segments(1).segment;
name = segments.segments(1).name;
[~, segment_m, indexes] = prepareInputForPlanner(segment);

%% single simulation
parameters.P_npDistances = [0.15, 0.2, 0.25];
parameters.P_ELDM = [-0.230023302976724	0.468296365511336	-0.190385591603131	-0.174774014076777	0.368721721474336	-0.211056276548872	-0.231339325522338;
-0.255003694198650	0.351441025614799	-0.0543410925120097	0.0271660325742101	-0.0792050028224978	0.0511467767037086	-0.231339325522338;
0.0216319451469262	-0.120626235909642	0.132344563327706	0.592550356037434	-0.870361141895546	0.278434887426610	-0.231339325522338];
%parameters.P_ELDM = reshape([2.86078399132926;0.884814215672794;-1.90657794718284;-3.09943416608130;-0.665457759838954;2.30236448840005;0.348462602099426;-0.107035325513227;-0.271014703397729;1.07959046302992;-0.775251579323662;-0.252977961446196;-0.822164501814478;1.36747233514778;0.113183483561418;-0.124241139196637;-0.454142531428492;0.293625990988783;-0.000983031283019174;-0.000983031283019174;-0.000983031283019174], 3,7);
%parameters.P_ELDM = zeros(3,7);

temp_folder_path = config.root;
plots_folder_name = 'plots';
set(0,'DefaultFigureVisible','off');

curveTypes = calculateCurveType(segment_m, indexes, name);

% oncoming traffic is scaled, with a low pass filter
oncomingTrafficScaled = segment_m(:,indexes.oncomingTraffic) > 0; 
oncomingTrafficScaled = movmean(oncomingTrafficScaled, 100); 
segment_m(:, end+1) = oncomingTrafficScaled;
indexes.oncomingTrafficScaled = size(segment_m,2);

% go on until stop line and find the optimal path
globalStartIndex = 2;
globalStopIndex = size(segment_m,1);
[path, U, dY, plannedPath, intentionPath, vehicleStateMemory] = pathGenerationLite();

end

function [path, U, dY, plannedPath, intentionPath, vehicleStateMemory] = pathGenerationLite()
    global segment_m indexes globalStartIndex globalStopIndex vehicleState parameters metadata x_k1 modelMode config
    % Entire model is calculated in the global frame!
    vehicleState = initVehicleState(segment_m, indexes, globalStartIndex, modelMode);
    vehicleStateCheck = initVehicleState(segment_m, indexes, globalStartIndex, "dynamicSimplified");
    replan = 1;
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
    initialized = false;

    try
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
                
                [coefficients, u, dy, X,Y,theta] = optimalPath (scenario, indexes, vehicleState, coefficients, c0, c1, parameters);
                if(initialized)
                    F = [F;getframe(gcf)];
                else
                    F = getframe(gcf);
                    initialized = true;
                end
                
    
                replan = 0;
                U = [U u]; dY = [dY dy];
                origo = [vehicleState.X, vehicleState.Y, vehicleState.theta]'; % in [m m rad]
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
            vehicleState = vehicleModel(vehicleState, dT, modelMode);            
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
    catch
        disp('Simulation failed on the way...');
    end

    writerObj = VideoWriter(fullfile(config.root, 'optimalPath.avi'));%% file name
    writerObj.FrameRate = 1/(0.05);
    fprintf('VIDEO GENERATION STARTED...\n');
    
    % open the video writer
    open(writerObj);
    % write the frames to the video
    for i=1:length(F)
        % convert the image to a frame
        frame = F(i) ;    
        writeVideo(writerObj, frame);
    end
    % close the writer object
    close(writerObj);
end

function [coefficients, u, dy, X,Y,theta, f] = optimalPath(scenario, indexes, vehicleState, coefficients, c0, c1, parameters)
% calculating cost sub-functions
% cost 1: instinctive path
[coefficients, u, dy, X,Y,theta] = eldmLite (scenario, indexes, vehicleState, coefficients, parameters.P_npDistances, parameters.P_ELDM, c0, c1);

% creating cost map from oncoming traffic
% it is defined by the furthermost np distance
% idea: the preview corridor is divided into grid members with the width of
% the corridor + thd with length of the furthermost node point
% first create the base map generated from the free driving path
dx = 0.05; % resolution in longitudinal direction is 5 cm
dy = 0.05; % resolution in lateral direction is 5 cm
ymax = 2; % 2 meters to left and right

% restoring waypoints
waypoints = [];
for i=1:length(coefficients)
    x = coefficients(i).sectionBorders(1):dx:coefficients(i).sectionBorders(2);
    waypoints = [waypoints; [x' (coefficients(i).coefficients(1) + coefficients(i).coefficients(2)*x + coefficients(i).coefficients(3)*x.^2 + coefficients(i).coefficients(4)*x.^3)']];
end
instinctivePathEgoFrame = waypoints;

% transforming corridor to the ego frame
x = 0:dx:coefficients(end).sectionBorders(2);
T = [cos(vehicleState.theta) sin(vehicleState.theta);
    -sin(vehicleState.theta) cos(vehicleState.theta)];
corridorGlobalFrame = [scenario(:,indexes.corrX) scenario(:,indexes.corrY)];
driverPathGlobalFrame = [scenario(:,indexes.X_abs) scenario(:,indexes.Y_abs)];
driverPathEgoFrame = (driverPathGlobalFrame - [vehicleState.X vehicleState.Y])*T';
corridorEgoFrame = (corridorGlobalFrame - [vehicleState.X vehicleState.Y])*T';
corridorEgoFrameResampled(:,2) = spline(corridorEgoFrame(:,1), corridorEgoFrame(:,2), x) + 1.875;
corridorEgoFrameResampled(:,3) = spline(corridorEgoFrame(:,1), corridorEgoFrame(:,2), x) - 1.875;
corridorEgoFrameResampled(:,1) = x;
centerLineEgoFrameResampled(:,2) = spline(corridorEgoFrame(:,1), corridorEgoFrame(:,2), x);
centerLineEgoFrameResampled(:,1) = x;

instinctivePathEgoFrameResampled(:,2) = spline(instinctivePathEgoFrame(:,1), instinctivePathEgoFrame(:,2), x);
instinctivePathEgoFrameResampled(:,1) = x;

instinctivePathEgoFrameResampledSpline(:,2) = spline(X,Y,x);
instinctivePathEgoFrameResampledSpline(:,1) = x;

nodePoints(:,1) = X; nodePoints(:,2) = Y;

% creating base map
a = 10; % second order coefficient for cost term
b = 0.1; % initial cost
for i=1:length(x)
    y_array(:,i) = corridorEgoFrameResampled(i,3):dy:corridorEgoFrameResampled(i,2); % from right lane edge to left lane edge
    y0 = instinctivePathEgoFrameResampled(i,2);
    costFreeDriving(:,i) = a*((y_array(:,i)-y0)/(2*ymax)).^2+b;
end

% getting the oncoming traffic information
c = 1; d = 1; k0 = 0.25; k1 = 2.25;
oncomingTrafficResampled = spline(corridorEgoFrame(:,1), scenario(:,indexes.oncomingTrafficScaled), x);
oncomingTrafficResampled = max(oncomingTrafficResampled,0);

for i=1:length(x)
    K = oncomingTrafficResampled(i); % between 0 and 1
    y = corridorEgoFrameResampled(i,2) - (corridorEgoFrameResampled(i,3):dy:corridorEgoFrameResampled(i,2)); 
    costTraffic(:,i) = (k0+K*k1)*exp(-(y-c))./(exp(-(y-c))+d);
end

% summing the cost map
w_fd = 0.5; w_tr = 0.5;
cost = w_fd*costFreeDriving + w_tr*costTraffic;

% plotting
f = figure();

h = pcolor(x, y_array, cost);
clim([0 1]); 
set(h,'LineStyle','none');
view(2);
xlim([0,coefficients(end).sectionBorders(2)]); ylim([min(-10, min(corridorEgoFrameResampled(:,3))),max(10, max(corridorEgoFrameResampled(:,2)))]);
hold on;
plot(instinctivePathEgoFrameResampled(:,1), instinctivePathEgoFrameResampled(:,2), 'LineWidth', 2, 'color', 'k');
plot(corridorEgoFrameResampled(:,1), corridorEgoFrameResampled(:,2), 'LineWidth', 2, 'LineStyle','--', 'color', 'g');
plot(corridorEgoFrameResampled(:,1), corridorEgoFrameResampled(:,3), 'LineWidth', 2, 'LineStyle','--', 'color', 'g');
plot(driverPathEgoFrame(:,1), driverPathEgoFrame(:,2), 'LineWidth', 1.5, 'color', 'r');
plot(centerLineEgoFrameResampled(:,1), centerLineEgoFrameResampled(:,2), 'LineWidth', 1.5, 'Color','w', 'LineStyle', '--');
plot(nodePoints(:,1), nodePoints(:,2), 'bo', 'LineWidth',3);
plot(instinctivePathEgoFrameResampledSpline(:,1), instinctivePathEgoFrameResampledSpline(:,2), 'LineWidth', 2, 'color', 'y');
%axis equal;
legend('cost map', 'instinctive path', 'lane left', 'lane right', 'driver path', 'centerline', 'node points', 'spline planned path');
colorbar; colormap(parula);
hold off;
pause(0.01);
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

    g = figure();

    plot(segment_m(:,indexes.X_abs), segment_m(:, indexes.Y_abs));
    grid on; hold on;
    plot(segment_m(curveTypes==1, indexes.X_abs), segment_m(curveTypes==1, indexes.Y_abs), 'ko');
    legend('original', 'straights');
    savefig(g, fullfile(temp_folder_path, plots_folder_name,...
                strcat('Map', name(1:end-4), '.fig')));
    saveas(g, fullfile(temp_folder_path, plots_folder_name,...
                strcat('Map_', name(1:end-4), '.png')));
    close(g);

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

    nominalSelect = "groundTruth";
    
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
end

function vehicleState = vehicleModel(vehicleState, dT, mode)
if(mode == "kinematic")    
    % kinematic bicycle model
    % parameters
    p_wheelBase = 2.7; % meter
    % calculating yaw-rate based on steering angle
    vehicleState.yawRate = tan(vehicleState.steeringAngle)/p_wheelBase * vehicleState.vx;
    vehicleState.theta = vehicleState.theta + dT*vehicleState.yawRate;
    % displacement
    dX = vehicleState.vx*dT*cos(vehicleState.theta);
    dY = vehicleState.vx*dT*sin(vehicleState.theta);
    % new vehicle position
    vehicleState.X = vehicleState.X+dX;
    vehicleState.Y = vehicleState.Y+dY;
    vehicleState.vy = 0;
    vehicleState.ax = 0;
    vehicleState.ay = vehicleState.vx*vehicleState.yawRate;
    
elseif (mode == "dynamic")    
    % dynamic bicycle model
    % only for high speeds
    % inputs:
    % - vehicleState.steeringAngle
    % - vehicleState.M_av/M_ah - breaking forces
    % - vehicleState.M_bv/M_bh - breaking forces
    % parameters:
    % r: radius of wheels in [m]
    % c_alfav/h: lateral slip coefficient [N/rad]
    % c_sv/h: longitudinal slip coefficient [N/%]
    % c_w: wind coefficient
    % rho_air: air density
    % A: preface
    % J: rotational intertia of the vehicle
    % m: mass of the vehicle
    % lv/lh: COG position from front and rear axle
    % Jwheel: rotiational interatia of the wheels
    r = 0.309725; c_alfaf = 76500; c_sf = 250; c_alfar = 76500; c_sr = 250;
    m = 1519;
    Jwheel = 250; J = 1818;
    A = 1.5; c_w = 0; rho_air = 1; lf = 1; lr = 1.5; 
    
%    vehicleState.s_v = -(vehicleState.v_vx_v - r*vehicleState.drho_v)/(max(abs(vehicleState.v_vx_v), r*vehicleState.drho_v)); % longitudinal slip front, wheel coordinate frame
    vehicleState.s_f = -(vehicleState.v_fx_v - r*vehicleState.drho_f)/(r*vehicleState.drho_f); % longitudinal slip front, wheel coordinate frame

%    vehicleState.alfa_v = -vehicleState.v_vy_v /(abs(r*vehicleState.drho_v)); % lateral slip front, wheel coordinate frame
    vehicleState.alfa_f = -vehicleState.v_fy_v /(r*vehicleState.drho_f); % lateral slip front, wheel coordinate frame

%    vehicleState.s_h = -(vehicleState.v_hx_v - r*vehicleState.drho_h)/(max(abs(vehicleState.v_hx_v), r*vehicleState.drho_h)); % longitudinal slip front, wheel coordinate frame
    vehicleState.s_r= -(vehicleState.v_rx_v - r*vehicleState.drho_r)/(r*vehicleState.drho_r); % longitudinal slip front, wheel coordinate frame

%    vehicleState.alfa_h = -vehicleState.v_hy_v /(abs(r*vehicleState.drho_h)); % lateral slip front, wheel coordinate frame
    vehicleState.alfa_r = -vehicleState.v_ry_v /(r*vehicleState.drho_r); % lateral slip front, wheel coordinate frame

    vehicleState.F_fx_v = c_sf*vehicleState.s_f; % longitudinal tyre force front, wheel coordinate frame
    vehicleState.F_fy_v = c_alfaf*vehicleState.alfa_f; % lateral tyre force front, wheel coordinate frame
    vehicleState.F_fx = cos(vehicleState.steeringAngle)*vehicleState.F_fx_v - sin(vehicleState.steeringAngle)*vehicleState.F_fy_v; % longitudinal tyre force front, vehicle frame
    vehicleState.F_fy = sin(vehicleState.steeringAngle)*vehicleState.F_fx_v + cos(vehicleState.steeringAngle)*vehicleState.F_fy_v; % lateral tyre force front, vehicle frame
    vehicleState.F_rx = c_sr*vehicleState.s_r; % longitudinal tyre force, rear wheel, vehicle coordinate frame == wheel coordinate frame
    vehicleState.F_ry = c_alfar*vehicleState.alfa_r; % lateral tyre force, rear wheel, vehicle coordinate frame == wheel coordinate frame
    vehicleState.F_wx = 0.5*c_w*rho_air*A*vehicleState.v_fx^2;
    vehicleState.F_wy = 0.5*c_w*rho_air*A*vehicleState.v_fy^2;
    
    % COG quantities    
    vehicleState.a_x = 1/m*(vehicleState.F_fx + vehicleState.F_rx - vehicleState.F_wx);
    vehicleState.a_y = 1/m*(vehicleState.F_fy + vehicleState.F_ry - vehicleState.F_wy);
    
    vehicleState.v_x = vehicleState.v_x + vehicleState.a_x*dT;
    vehicleState.v_y = vehicleState.v_y + vehicleState.a_y*dT;
    
    vehicleState.eps_sigma = 1/J*(lf*vehicleState.F_fy - lr*vehicleState.F_ry);
    vehicleState.yawRate = vehicleState.yawRate + dT*vehicleState.eps_sigma;
    
    % Baselink quantities
    vehicleState.ax = vehicleState.a_x;
    vehicleState.vx = vehicleState.v_x;
    vehicleState.ay = vehicleState.yawRate*vehicleState.vx;
    vehicleState.vy = 0;
    
    % transforming to front and rear wheel
    vehicleState.v_fx = vehicleState.v_x;
    vehicleState.v_fy = vehicleState.v_y + vehicleState.yawRate*lf;
    vehicleState.v_rx = vehicleState.v_x;
    vehicleState.v_ry = vehicleState.v_y - vehicleState.yawRate*lr;
    
    vehicleState.v_fx_v = cos(vehicleState.steeringAngle)*vehicleState.v_fx + sin(vehicleState.steeringAngle)*vehicleState.v_fy;
    vehicleState.v_fy_v = -sin(vehicleState.steeringAngle)*vehicleState.v_fx + cos(vehicleState.steeringAngle)*vehicleState.v_fy;
    vehicleState.v_rx_v = vehicleState.v_rx;
    vehicleState.v_ry_v = vehicleState.v_ry;
    
    vehicleState.ddrho_f = 1/Jwheel * (vehicleState.M_af - vehicleState.M_bf)*sign(vehicleState.drho_f) - r*vehicleState.F_fx_v;
    vehicleState.ddrho_r = 1/Jwheel * (vehicleState.M_ar - vehicleState.M_br)*sign(vehicleState.drho_r) - r*vehicleState.F_rx;
    
    vehicleState.drho_f = vehicleState.drho_f + vehicleState.ddrho_f*dT;
    vehicleState.drho_r = vehicleState.drho_r + vehicleState.ddrho_r*dT;
    
    % absolute frame quantities
    vehicleState.theta = vehicleState.theta + dT*vehicleState.yawRate;
    % displacement
    dX = vehicleState.v_rx*dT*cos(vehicleState.theta)-sin(vehicleState.theta)*vehicleState.v_ry*dT;
    dY = vehicleState.v_rx*dT*sin(vehicleState.theta)+cos(vehicleState.theta)*vehicleState.v_ry*dT;
    % new vehicle position
    vehicleState.X = vehicleState.X+dX;
    vehicleState.Y = vehicleState.Y+dY;
    
    vehicleState.v = (vehicleState.v_x^2 + vehicleState.v_y^2)^0.5;  
    
elseif (mode == "dynamicSimplified")
    r = 0.309725; c_alfaf = 76500; c_sf = 250; c_alfar = 76500; c_sr = 250;
    m = 1519;
    Jwheel = 250; J = 1818;
    lf = 1; lr = 1.5;     
    
    df = vehicleState.steeringAngle;
    vehicleState.F_fx = cos(df)*(-c_sf*(cos(df)*vehicleState.v_x+sin(df)*vehicleState.v_y+sin(df)*vehicleState.yawRate*lf-r*vehicleState.drho_f)/(r*vehicleState.drho_f))+ ...
        sin(df)*c_alfaf*(-sin(df)*vehicleState.v_x+cos(df)*vehicleState.v_y+cos(df)*vehicleState.yawRate*lf)/(r*vehicleState.drho_f);
    vehicleState.F_fy = -sin(df)*c_sf*(cos(df)*vehicleState.v_x+sin(df)*vehicleState.v_y+sin(df)*vehicleState.yawRate*lf-r*vehicleState.drho_f)/(r*vehicleState.drho_f) - ...
        cos(df)*c_alfaf*(-sin(df)*vehicleState.v_x+cos(df)*vehicleState.v_y+cos(df)*vehicleState.yawRate*lf)/(r*vehicleState.drho_f);
    vehicleState.F_rx = c_sr*(-(vehicleState.v_x-r*vehicleState.drho_r)/(r*vehicleState.drho_f));
    vehicleState.F_ry = c_alfar*(-(vehicleState.v_y-vehicleState.yawRate*lr)/(r*vehicleState.drho_r));
    vehicleState.a_x = 1/m*(vehicleState.F_fx + vehicleState.F_rx);
    vehicleState.a_y = 1/m*(vehicleState.F_fy + vehicleState.F_ry);
    vehicleState.v_x = vehicleState.v_x + vehicleState.a_x*dT;
    vehicleState.v_y = vehicleState.v_y + vehicleState.a_y*dT;
    vehicleState.yawRate = vehicleState.yawRate + dT/J*lf*vehicleState.F_fy - dT*lr/J*vehicleState.F_ry;
    vehicleState.drho_f = vehicleState.drho_f+dT*r*c_sf*(cos(df)*vehicleState.v_x+sin(df)*vehicleState.v_y+sin(df)*vehicleState.yawRate*lf-r*vehicleState.drho_f)/(r*vehicleState.drho_f);
    vehicleState.drho_r = vehicleState.drho_r+r*dT*c_sr*(vehicleState.vx-r*vehicleState.drho_r)/(r*vehicleState.drho_r);
    
    vehicleState.theta = vehicleState.theta + vehicleState.yawRate*dT;
    % displacement
    v_hx = vehicleState.v_x;
    v_hy = vehicleState.v_y - lr*vehicleState.yawRate;
    dX = v_hx*dT*cos(vehicleState.theta)-sin(vehicleState.theta)*v_hy*dT;
    dY = v_hx*dT*sin(vehicleState.theta)+cos(vehicleState.theta)*v_hy*dT;
    % new vehicle position
    vehicleState.X = vehicleState.X+dX;
    vehicleState.Y = vehicleState.Y+dY;
    
end
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
