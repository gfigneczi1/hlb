clear; close all;

global modelID path GENERATE_OPT_PLOTS globalStartIndex globalStopIndex vehicleState metadata  modelID modelMode
modelID = "eldm";

parameters.P_npDistances = [0.04, 0.116, 0.388];
parameters.numberOfNodePoints = 3;
parameters.drID = 11;
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

GENERATE_MAP_PLOTS = false;
GENERATE_OPT_PLOTS = false;

Ts = 0.02;
T = 30;
N = T/Ts;

% generate parameters
driverData = load("C:\git\KDP\HLB_for_AV_MotionControl\02_Results\MassMeasurements\summaryData\Dr023_2023-05-22_12-17-22_r31_withTraffic.mat");
fn = fieldnames(driverData);
for i=1:length(fn)
    driverData.(fn{i}) = driverData.(fn{i})';
end

[~, segment_m, indexes] = prepareInputForPlanner(struct2table(driverData));
             
globalStartIndex = 2; % minimum is 2, otherwise it fails!
globalStopIndex = size(segment_m,1);
             
% Definition of metadata
metadata.pathValidity = 1;
modelID = "groundTruth";
modelMode = "kinematic";
        
[~, U, dY, ~, intentionPath] = functional_pathGeneration(segment_m, indexes, parameters);
P_ELDM = functional_driverModelLearning(U', dY' ,8);
% parameters.P_ELDM = P_ELDM;
modelID = "eldm";
        
%[~, U, dY, ~, intentionPath] = functional_pathGeneration(segment_m, indexes, parameters);

path = load("C:\git\hlb\Solutions\CruiseConcept\modelPredictiveControl\31_south_map.mat");
path = path.traj_local;
path = path(200:end,:);
orientation = diff(path(:,2))./diff(path(:,1));
orientation = [orientation; orientation(end)];
orientation = atan(movmean(orientation,25));
[poses,directions] = smoothPathSpline([path orientation],ones(size(path,1),1),floor(size(path,1)/4));
path = [path(:,1), spline(poses(:,1),poses(:,2),path(:,1))];
path(:,3) = orientation;
path(:,4) = zeros(size(path,1),1);
indexes.LaneOrientation = 4;
curvature = diff(orientation)./diff(path(:,1));
curvature = [curvature; curvature(end)];
curvature = movmean(curvature,50);
path(:,5) = curvature;
indexes.X_abs = 1;
indexes.Y_abs = 2;
indexes.theta_calc = 3;
indexes.corrX = 1;
indexes.corrY = 2;
indexes.LaneCurvature = 5;

global x0 parameters pointsRefGlobal pointsOppositeLaneGlobal orientationRefVectorGlobal trafficDistance map priorPath partialCosts ay0
parameters.Tsolve = 0.02; % only relevant when self-implemented solver is used
parameters.vx = 17; % m/s

parameters.Lr = 1.8; % m
parameters.Lf = 1.2; % m

parameters.alfa = [0.25 0; 0 0.0];
parameters.kappa = [1e-2, 0.0];
parameters.beta = 0.0; %0.1
parameters.Sresidual = 0.0;
parameters.delta = 0.0; %0.8
parameters.maximumRelevantTrafficTime = 4;

parameters.Th = 1.5;
parameters.Np = 20;
parameters.aymax = 2.5;
parameters.umax = 5*pi()/180;
parameters.jymax = 1;
parameters.allowedDeviation = 0.0;
parameters.laneWidth = 4; %m

parameters.costMode = "costMap";

% planner
parameters.P_npDistances = [0.04, 0.116, 0.388];
parameters.P_ELDM = [0.39 0.24 0.12 -0.23 -0.39 -0.18 0.07; ...
    -0.39 -0.2 -0.1 0.43 0.6 0.18 0.07; ...
    0.08 0.04 0.07 -0.18 -0.18 0.03 0.07];
%parameters.P_ELDM = [reshape(P_ELDM(1:9),3,3) reshape(P_ELDM(10:18),3,3) ones(3,1)*P_ELDM(19)];
%parameters.P_ELDM = zeros(3,7);

parameters.traffic = generateTraffic(parameters.vx, 50, path, T, Ts);
    
U = zeros(1,parameters.Np);
u_saved = [];

Td = parameters.Th/parameters.Np;

M = 0;

x_k = zeros(3,1);
x_k(1) = path(1,1); x_k(2) = path(1,2); x_k(3) = path(1,3);
vehicleState.X = x_k(1);
vehicleState.Y = x_k(2);
vehicleState.theta = x_k(3);

egoPath = [vehicleState.X vehicleState.Y];
globalPosteriorPath = egoPath;
% simulation loop
try
for k=1:N
    
    x0 = x_k;
    t_k = k*Ts;
    if (true) %(t_k/(M*Td) >= 1)
        % new plan - store previous global map points up until now
        if (k>1)
            globalPosteriorPath = [globalPosteriorPath; [vehicleState.X spline(priorPath(:,1), priorPath(:,2), vehicleState.X)]];
        end

        % generate local map
        nearestIndex = getNearestIndex(path(:,1),path(:,2), [vehicleState.X vehicleState.Y]);
        scenario = functional_cutScenario(path, indexes, nearestIndex, vehicleState, [0,160]);
        priorPath = functional_planner (scenario, indexes, globalPosteriorPath, parameters);
        % finding local path, interpolate to prediction horizon
        distances = (sum((priorPath(:,1:2)-[x_k(1) x_k(2)])'.^2)).^0.5;
        nearestIndex = find(distances==min(distances),1);
        if (~isempty(nearestIndex))
            % there is a valid point found globally, transform trajectory to
            % local
            dx = parameters.vx*parameters.Th/parameters.Np;
            xend = dx*parameters.Np; % distance during one optimization cycle
            T = [cos(x_k(3)) sin(x_k(3)); -sin(x_k(3)) cos(x_k(3))];
            P = [x_k(1) x_k(2)];
            localPath = (priorPath(nearestIndex:end,1:2)-P)*T';
            localDistances = (sum(localPath'.^2)).^0.5;
            endPoint = find(localDistances(1:end)>=xend,1);
            xRefVector = linspace(dx,xend, parameters.Np);
            yRefVector = spline(localPath(1:endPoint,1), localPath(1:endPoint,2), xRefVector);
            orientationRefVector = movmean(diff(yRefVector)./diff(xRefVector), 10);
            orientationRefVector = [orientationRefVector orientationRefVector(end)];
            % transform back to global
            pointsRefGlobal = [xRefVector' yRefVector']*T+P;
            orientationRefVectorGlobal = orientationRefVector'+x_k(3);

            % parameters.traffic have the same time scaled traffic information
            % cut the horizon from the traffic information
            t00 = find(parameters.traffic(:,end)>=t_k,1);
            t01 = find(parameters.traffic(:,end)>=t_k+parameters.Th,1);
            traffic = parameters.traffic(t00:t01,:);
            % resample to the horizon length
            timeResampled = t_k+(Td:Td:parameters.Th);
            for i=1:size(traffic,2)-1
                trafficResampled(:,i) = spline(traffic(:,end), traffic(:,i), timeResampled');
            end
            trafficResampled(:,end+1) = timeResampled';
            % calculate traffic distance
            trafficDistance = ((trafficResampled(:,1)-pointsRefGlobal(:,1)).^2+(trafficResampled(:,2)-pointsRefGlobal(:,2)).^2).^0.5;
            dlim = parameters.maximumRelevantTrafficTime*parameters.vx*2;
            trafficDistance = trafficDistance/dlim;
        end
        % finding local path, interpolate to prediction horizon
        trafficReordered = parameters.traffic(end:-1:1,1:2);
        distances = (sum((trafficReordered-[x_k(1) x_k(2)])'.^2)).^0.5;
        nearestIndex = find(distances==min(distances),1);
        if (~isempty(nearestIndex))
            % there is a valid point found globally, transform trajectory to
            % local
            dx = parameters.vx*parameters.Th/parameters.Np;
            xend = dx*parameters.Np; % distance during one optimization cycle
            T = [cos(x_k(3)) sin(x_k(3)); -sin(x_k(3)) cos(x_k(3))];
            P = [x_k(1) x_k(2)];
            localPath = (trafficReordered(nearestIndex:end,1:2)-P)*T';
            localDistances = (sum(localPath'.^2)).^0.5;
            endPoint = find(localDistances(1:end)>=xend,1);
            if(isempty(endPoint) || size(localPath,1)<25)
                parameters.delta = 0;
                yRefVector = xRefVector*0;
            else
                yRefVector = spline(localPath(1:endPoint,1), localPath(1:endPoint,2), xRefVector);
            end
            % transform back to global
            pointsOppositeLaneGlobal = [xRefVector' yRefVector']*T+P;
        end
        %U = zeros(1,parameters.Np);
        M = M+1;
        if(~isempty(u_saved))
            ay0 = parameters.vx^2*tan(u_saved(end))/(parameters.Lr+parameters.Lf);
        else
            ay0 = 0;
        end
        for o=1:1
            [U, fval(M,1), exitFlag] = fminunc(@objfunx2, U); 
        end
        fval(M,2:6) = partialCosts;
        for i=1:parameters.Np
            t_i = i*Td;
            u_i = mu(U, t_i);
            if (i==1)
                x(:,i) = kinematicSingleTrack(u_i, x_k, Td);
            else
                x(:,i) = kinematicSingleTrack(u_i, x(:,i-1), Td);
            end  
        end 
    end
    u_k = mu(U, t_k-(M-1)*Td); % receeding horizon

    % saving for debugging
    u_saved(k) = u_k;
    x_saved(k,1:size(x_k,1)) = x_k';
    th_saved(k) = t_k-(M-1)*Td;
    kappa_saved(k) = scenario(1,indexes.LaneCurvature);

    % state update
    x_k = kinematicSingleTrack(u_k, x_k, Ts);
    T = [cos(x_k(3)) sin(x_k(3)); -sin(x_k(3)) cos(x_k(3))];
    localPath = (path(:,1:2)-[vehicleState.X vehicleState.Y])*T';

    error(k) = -spline(localPath(:,1), localPath(:,2), 0);

    localPriorPath = (priorPath(:,1:2)-[vehicleState.X vehicleState.Y])*T';
    errorToInstinctivePath(k) = -spline(localPriorPath(:,1), localPriorPath(:,2), 0);
    errorInstinctivePathToPath(k) = error(k)-errorToInstinctivePath(k);
    
    if(GENERATE_MAP_PLOTS)
        figure(2);
        subplot(3,1,1);
        plot(x_saved(k,1), x_saved(k,2), 'kx', x(1, :), x(2, :), 'r--', priorPath(:,1), priorPath(:,2), 'b:', globalPosteriorPath(:,1), globalPosteriorPath(:,2), 'ko'); hold on; grid on;
        title(strcat("error is:", num2str(error(k))));
        subplot(3,1,2);
        plot(x_saved(k,1), u_saved(k), 'kx', x(1, :), U, 'r--'); hold on; grid on; xlabel('X(m)'); ylabel('\delta_f')
        subplot(3,1,3);
        plot(x_saved(k,1), error(k), 'kx'); grid on; hold on; ylabel('error(m)');
        pause(0.1);
    end

    % update vehicle state struct
    vehicleState.X = x_k(1);
    vehicleState.Y = x_k(2);
    vehicleState.theta = x_k(3);

    egoPath = [egoPath; [vehicleState.X vehicleState.Y]];
    clc;
    fprintf("Iteration number %d / %d and the error is %.3f and error statistics: %.3f, %.3f, %.3f\n", k, N, errorToInstinctivePath(k), max(abs(errorToInstinctivePath)),  mean(abs(errorToInstinctivePath)),std(errorToInstinctivePath));
    fprintf("\tTotal cost=%f,\n\toutput_cost=%f,\n\tinput cost=%f,\n\tinput gradient cost is %f\n\ttraffic cost: %f\n\trun off cost: %f\n", fval(M,1), fval(M,2), fval(M,3), fval(M,4), fval(M,5), fval(M,6));
    fprintf("Input normalized is %f\n", u_k(1)/parameters.umax);
end
catch e
    N = length(error);
end

f = figure();
set(f,'defaulttextInterpreter','latex') ;
set(f, 'defaultAxesTickLabelInterpreter','latex');  
set(f, 'defaultLegendInterpreter','latex');
set(f, "Position", [100,100,505,750]);
subplot(4,1,1);
% error plot
plot(Ts:Ts:N*Ts, errorInstinctivePathToPath, 'DisplayName','Planned lane offset', 'LineWidth',1.5, 'Color','k', 'LineStyle','--');
grid on; hold on;
plot(Ts:Ts:N*Ts, error, 'DisplayName','Real lane offset', 'LineWidth',1.5,'Color','k');
ylim([-1,1]); xlabel('time(s)'); ylabel('$\delta(t)$');
xline(0,'Alpha',0.8,'Color','k', 'HandleVisibility','off');
xlim([Ts,N*Ts]);
set(gca,'FontSize',14);
legend('Location','best', 'FontSize',11);
if (parameters.costMode=="normal" && parameters.delta==0 && all(all(parameters.P_ELDM==0)))
    title("\textbf{Normal mode, centerline driving}");
elseif (parameters.costMode=="costMap" && parameters.delta==0 && ~all(all(parameters.P_ELDM==0)))
    title("\textbf{Modified mode, instinctive path}");
elseif (parameters.costMode=="costMap" && parameters.delta>0 && ~all(all(parameters.P_ELDM==0)))
    title("\textbf{Modified mode, instinctive+reactive path}");
end

subplot(4,1,2);
% input plot
plot(Ts:Ts:N*Ts, u_saved, 'DisplayName', 'Steering angle', 'LineWidth',1.5, 'color', 'k');
grid on;
xlabel('time(s)'); ylabel('$\alpha_f(rad)$');
ylim([-0.05,0.05]); xlim([Ts,N*Ts]);
xline(0,'Alpha',0.8,'Color','k', 'HandleVisibility','off');
set(gca,'FontSize',14);
legend('Location','best', 'FontSize',11);

subplot(4,1,3);
% cost plot
plot(Ts:Ts:N*Ts, fval(:,1), 'DisplayName', 'Total cost', 'LineWidth',2,'color','k');
hold on; grid on;
plot(Ts:Ts:N*Ts, fval(:,2), 'DisplayName', 'Output error cost', 'LineWidth',1.5,'color','k', 'LineStyle', '--');
plot(Ts:Ts:N*Ts, fval(:,3), 'DisplayName', 'Input cost', 'LineWidth',1.5,'color','b', 'LineStyle','--', 'LineStyle', ':');
plot(Ts:Ts:N*Ts, fval(:,4), 'DisplayName', 'Input gradient cost', 'LineWidth',1.5,'color','k','LineStyle', '-.');
plot(Ts:Ts:N*Ts, fval(:,5), 'DisplayName', 'Traffic cost', 'LineWidth',1.5,'color','b','LineStyle', '-');
plot(Ts:Ts:N*Ts, fval(:,6), 'DisplayName', 'Run-off cost', 'LineWidth',1.5,'color','r', 'LineStyle', '-');
ylim([0,1]);
xlabel('time(s)'); ylabel('cost(-)');
xlim([Ts,N*Ts]);
set(gca,'FontSize',14);
legend('Location','north', 'FontSize',11, 'Orientation','horizontal','NumColumns',2);

subplot(4,1,4);
% curvature plot
plot(Ts:Ts:N*Ts, kappa_saved, 'DisplayName', 'Road curvature', 'Color','k', 'LineWidth',1.5);
hold on; grid on;
xline(0,'Alpha',0.8,'Color','k', 'HandleVisibility','off');
xlabel('time(s)'); ylabel('$\kappa_{road}(1/m)$');
ylim([-8e-3,8e-3]);
xlim([Ts,N*Ts]);
set(gca,'FontSize',14);
legend('Location','best', 'FontSize',11);


function f = objfunx2(U)
global x0 parameters pointsRefGlobal pointsOppositeLaneGlobal trafficDistance map priorPath path orientationRefVectorGlobal GENERATE_OPT_PLOTS partialCosts ay0 t0
    
    Td = parameters.Th/parameters.Np;
    x_k = x0;
    for k=1:parameters.Np
        t_k = k*Td;
        u_k = mu(U, t_k);
        x_k = kinematicSingleTrack(u_k, x_k, Td);
        x(:, k) = x_k;
    end 
    % transform everything to local map
    T = [cos(x0(3)) sin(x0(3)); -sin(x0(3)) cos(x0(3))];
    xLocal(:,1:2) = (x(1:2,:)'-x0(1:2)')*T';
    xLocal(:,3) = x(3,:)'-x0(3);
    pointsRefLocal = (pointsRefGlobal-x0(1:2)')*T';
    orientationRefVectorLocal = orientationRefVectorGlobal-x0(3);
    pointsOppositeLaneLocal = (pointsOppositeLaneGlobal-x0(1:2)')*T';
    % output
    y = [0 1 0; 0 0 1]*xLocal';
    % cost 
    if(parameters.costMode=="costMap")
        f = 0;
        yref = [pointsRefLocal(:,2)'; orientationRefVectorLocal'];
        for i=1:size(x,2)
            if(abs(y(2,i)-yref(1,i))<parameters.allowedDeviation)
                df=0;
            else
                dist = y(1,i)-yref(1,i);
                df = parameters.alfa(1,1)*(kernelInverseGaussian(dist/0.5, 0.5, 0) - ...
                    kernelInverseGaussian(0, 0.5, 0));
                %df = parameters.alfa(1,1)*exp(0.25*abs(y(1,i)-yref(1,i))/(0.5*parameters.laneWidth));
            end
            f = f+df;
        end
        f = f/size(y,2);
        partialCosts(1) = f;
        % input error
        % recalculation to lateral acceleration - ackermann angle
        Ay = parameters.vx^2*tan(U)/(parameters.Lr+parameters.Lf);
        f = f+(Ay/parameters.aymax*parameters.kappa(1)*Ay'/parameters.aymax)/size(y,2);
        partialCosts(2) = (Ay/parameters.aymax*parameters.kappa(1)*Ay'/parameters.aymax)/size(y,2);
        % input gradient error
        Aymod = [ay0 Ay];
        dAy = diff(Aymod);
        Jy = dAy/Td;
        f = f+(Jy/parameters.jymax*parameters.kappa(2)*Jy'/parameters.jymax)/size(y,2);
        partialCosts(3) = (Jy/parameters.jymax*parameters.kappa(2)*Jy'/parameters.jymax)/size(y,2);
    
        % traffic cost
        if (min(trafficDistance) < 1.5)
            costTraffic = kernelGaussian(trafficDistance, 0.75, 0);
            fTraffic = 0;
            for i=1:size(x,2)
                distToOpponentLane = pointsOppositeLaneLocal(i,2)-y(1,i);
                df = costTraffic(i)*(1-sigmoid(distToOpponentLane/(parameters.laneWidth*1.5),6));
                fTraffic = fTraffic+df;
            end
            % calculating the distance to the opposite lane centerline
            f = f+parameters.delta*fTraffic;
            partialCosts(4) = parameters.delta*fTraffic;
        else
            partialCosts(4) = 0;
        end

        % right hand side run-off cost
        fRunOff = 0;
        for i=1:size(x,2)
            distToRightSide = parameters.laneWidth/2-(pointsOppositeLaneLocal(i,2)-parameters.laneWidth-y(1,i));
            df = (1-sigmoid(distToRightSide/(parameters.laneWidth*0.5),10));
            fRunOff = fRunOff+df;
        end
        f = f+parameters.beta*fRunOff;
        partialCosts(5) = parameters.beta*fRunOff;
    else
        yref = [yRefVectorGlobal'; orientationRefVectorGlobal'];
        f = 0;
        % accumulated error
        for n=1:size(y,1)
            yn = y(n,:); yrefn = yref(n,:);
            if (n==1)
                normFactor = 0.5;
            else
                normFactor = pi()/12;
            end
            f = f+(yn-yrefn)/normFactor*parameters.alfa(n,n)*(yn-yrefn)'/normFactor;
        end
        f = f/size(y,2);
        partialCosts(1) = f;
        % input error
        % recalculation to lateral acceleration - ackermann angle
        Ay = parameters.vx^2*tan(U)/(parameters.Lr+parameters.Lf);
        f = f+(Ay/parameters.aymax*parameters.kappa(1)*Ay'/parameters.aymax)/size(y,2);
        partialCosts(2) = (Ay/parameters.aymax*parameters.kappa(1)*Ay'/parameters.aymax)/size(y,2);
        partialCosts(3:4) = nan;
    end
    
    if (GENERATE_OPT_PLOTS)
        figure(1)
        plot(x(1,:), x(2,:), 'r', x(1,:), yRefVectorGlobal(1:end), 'k', priorPath(:,1), priorPath(:,2), 'g', path(:,1), path(:,2), 'k:', 'LineWidth',2);
        title(strcat("Total cost=", num2str(f), ",outputCost=", num2str(partialCosts(1)), ",inputCost=", num2str(partialCosts(2))));
        xlim([min(x(1,:))-3, max(x(1,:))+3]);
        ylim([min(x(2,:))-1, max(x(2,:))+1]);
        pause(0.05);
    end

end

function u_k = mu(U, t_k)
    global parameters
    Td = parameters.Th/parameters.Np;
    N = max(1,floor(t_k/Td));
    u_k = U(N);
end

function cost = kernelInverseGaussian(x, sigma, mu)
    a = 1/(sigma*sqrt(2*pi()));
    b = mu;
    c = sigma;
    cost = a*exp(-(x-b).^2/(2*c^2))/(a*exp(-(0).^2/(2*c^2)));
    cost = 1-cost;
    %cost = cost/(1-min(cost));
    %cost = cost-min(cost); % shift to zero local extremum    
end

function cost = kernelGaussian(x, sigma, mu)
    a = 1/(sigma*sqrt(2*pi()));
    b = mu;
    c = sigma;
    cost = a*exp(-(x-b).^2/(2*c^2))/(a*exp(-(0).^2/(2*c^2)));
end

function x_k = kinematicSingleTrack(u_k, x_k1, Ts)
global parameters
x_k(1,1) = x_k1(1,1)+Ts*parameters.vx*cos(x_k1(3));
x_k(2,1) = x_k1(2,1)+Ts*parameters.vx*sin(x_k1(3));
x_k(3,1) = x_k1(3,1)+parameters.vx*Ts*tan(u_k(1))/(parameters.Lr+parameters.Lf);
end

function y = sigmoid(x,a)
    % input x is normalized to 0 and 1, therefore its multiplication with a
    % is needed
    xmax = 16/a;
    xn = x*xmax;
    y= 1./(1+exp(-a*(xn-1)));
end

function traffic = generateTraffic(speed, offset, path, Th, Ts)
global parameters
% speed: constant speed of contraflow traffic
% offset: in time, compared to the other end of the road
% contraflow traffic two, with the help of the ELDM model structure
% path: path of the ego vehicle
% T: length of measurement (s), sample time of measurement (s)
positionOffset = offset*speed; % in meters

for i=1:size(path,1)
    T = [cos(path(i,3)) -sin(path(i,3)); sin(path(i,3)) cos(path(i,3))];
    oppositeLanePath(i,1:2) = path(i,1:2)+[0 parameters.laneWidth]*T';
end
oppositeLanePath = oppositeLanePath(end:-1:1,:);
orientation = diff(oppositeLanePath(:,2))./diff(oppositeLanePath(:,1));
orientation = [orientation; orientation(end)];
orientation = atan(movmean(orientation,25))+pi();
oppositeLanePath(:,3) = orientation;
oppositeLanePath(:,4) = zeros(size(oppositeLanePath,1),1);
curvature = diff(oppositeLanePath(:,3))./abs(diff(oppositeLanePath(:,1)));
curvature = [curvature; curvature(end)];
curvature = movmean(curvature,50);
oppositeLanePath(:,5) = curvature;
indexes.X_abs = 1;
indexes.Y_abs = 2;
indexes.theta_calc = 3;
indexes.LaneOrientation = 4;
indexes.corrX = 1;
indexes.corrY = 2;
indexes.LaneCurvature = 5;

% cut data from position offset
runDistances = cumtrapz((diff(oppositeLanePath(:,1)).^2+diff(oppositeLanePath(:,2)).^2).^0.5);
idx = find(runDistances>positionOffset,1);
oppositeLanePath = oppositeLanePath(idx:end,:);

% resample due to speed
p0 = oppositeLanePath(1,:);
time = 0:Ts:Th;
dx = speed*Ts;
oppositeLanePathResampled(1,:) = p0;
for i=2:length(time)
    p1 = [p0(1)+dx*cos(p0(3)) p0(2)+dx*sin(p0(3))];
    oppositeLanePathResampled(i,1) = p1(1);
    oppositeLanePathResampled(i,2) = spline(oppositeLanePath(:,1), oppositeLanePath(:,2), p1(1));
    oppositeLanePathResampled(i,3) = spline(oppositeLanePath(:,1), oppositeLanePath(:,3), p1(1));
    oppositeLanePathResampled(i,4) = spline(oppositeLanePath(:,1), oppositeLanePath(:,4), p1(1));
    oppositeLanePathResampled(i,4) = spline(oppositeLanePath(:,1), oppositeLanePath(:,5), p1(1));
    p1(2) = oppositeLanePathResampled(i,2);
    p1(3) = oppositeLanePathResampled(i,3);
    p0 = p1;
end

oppositeLanePathResampled(:,end+1) = time(1:end)';

traffic = oppositeLanePathResampled;
end