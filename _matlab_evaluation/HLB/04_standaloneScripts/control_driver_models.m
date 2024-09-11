clear; close all;

global modelID path GENERATE_OPT_PLOTS globalStartIndex globalStopIndex vehicleState metadata  modelID modelMode parameters
modelType = "mpc";

parameters.drID = 11;
parameters.vehicleParameters.wheelBase = 2.7;
parameters.vehicleParameters.r = 0.309725; 

parameters.pp.lat = 2.25; % seconds
parameters.stanley.k = 0.5;
parameters.pid.P = 0.01;
parameters.pid.I = 0.01;
parameters.pid.D = 0.05;

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
                         
% Definition of metadata
modelMode = "kinematic";

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

global x0 pointsRefGlobal pointsOppositeLaneGlobal yRefVector orientationRefVector orientationRefVectorGlobal trafficDistance map priorPath partialCosts ay0
parameters.Tsolve = 0.02; % only relevant when self-implemented solver is used
parameters.vx = 17; % m/s

parameters.Lr = 1.8; % m
parameters.Lf = 1.2; % m

parameters.alfa = [0.0 0; 0 2];
parameters.kappa = [2, 0.0];
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
int_e = 0;
e_k1 = 0;
% simulation loop
for k=1:N
    
    x0 = x_k;
    t_k = k*Ts;
    % finding local path, interpolate to prediction horizon
    distances = (sum((path(:,1:2)-[x_k(1) x_k(2)])'.^2)).^0.5;
    nearestIndex = find(distances==min(distances),1);
    if (~isempty(nearestIndex))
        % there is a valid point found globally, transform trajectory to
        % local
        dx = parameters.vx*parameters.Th/parameters.Np;
        xend = dx*parameters.Np; % distance during one optimization cycle
        T = [cos(x_k(3)) sin(x_k(3)); -sin(x_k(3)) cos(x_k(3))];
        P = [x_k(1) x_k(2)];
        localPath = (path(nearestIndex:end,1:2)-P)*T';
        localDistances = (sum(localPath'.^2)).^0.5;
        endPoint = find(localDistances(1:end)>=xend,1);
        xRefVector = linspace(dx,xend, parameters.Np);
        yRefVector = spline(localPath(1:endPoint,1), localPath(1:endPoint,2), xRefVector);
        orientationRefVector = movmean(diff(yRefVector)./diff(xRefVector), 10);
        orientationRefVector = [orientationRefVector orientationRefVector(end)];
        % transform back to global
        pointsRefGlobal = [xRefVector' yRefVector']*T+P;
        orientationRefVectorGlobal = orientationRefVector'+x_k(3);
    end
    M = M+1;
    if(~isempty(u_saved))
        ay0 = parameters.vx^2*tan(u_saved(end))/(parameters.Lr+parameters.Lf);
    else
        ay0 = 0;
    end
    tic;
    switch modelType
        case "mpc"
            for o=1:1
                [U, fval(M,1), exitFlag] = fminunc(@objfunx2, U); 
            end
            for i=1:parameters.Np
                t_i = i*Td;
                u_i = mu(U, t_i);
                if (i==1)
                    x(:,i) = kinematicSingleTrack(u_i, x_k, Td);
                else
                    x(:,i) = kinematicSingleTrack(u_i, x(:,i-1), Td);
                end  
            end
            u_k = mu(U, t_k-(M-1)*Td); % receeding horizon
        case "pp"
            % calculate the preview point
            y = spline(localPath(:,1), localPath(:,2), parameters.pp.lat*parameters.vx);
            L = sqrt((parameters.pp.lat*parameters.vx)^2+y^2);
            r = L^2/(2*y);
            u_k = atan((parameters.Lr+parameters.Lf)/r);
        case "stanley"
            % cross track error
            p0 = [0 spline(localPath(:,1), localPath(:,2), 0)];
            p1 = [(parameters.Lr+parameters.Lf) spline(localPath(:,1), localPath(:,2), (parameters.Lr+parameters.Lf))];
            e = p1(2);
            % cross track steering
            u_k_ct = atan(parameters.stanley.k*e/parameters.vx);
            % heading error
            theta = atan((p1(2)-p0(2))/(p1(1)-p0(1)));
            % total steering input
            u_k = u_k_ct+theta;
        case "pid"
            p1 = [0 spline(localPath(:,1), localPath(:,2), 0)];
            e = p1(2);
            int_e = int_e+e*Ts;
            der_e = (e-e_k1)/Ts;
            u_k = parameters.pid.P*e+parameters.pid.I*int_e+parameters.pid.D*der_e;
            e_k1 = e;
    end
    t = toc;

    % state update
    x_k = kinematicSingleTrack(u_k, x_k, Ts);
    T = [cos(x_k(3)) sin(x_k(3)); -sin(x_k(3)) cos(x_k(3))];
    localPath = (path(:,1:2)-[vehicleState.X vehicleState.Y])*T';

    error(k) = -spline(localPath(:,1), localPath(:,2), 0);

    % saving for debugging
    u_saved(k) = u_k;
    x_saved(k,1:size(x_k,1)) = x_k';
    th_saved(k) = t_k-(M-1)*Td;
    compt_saved(k) = t;
    
    if(GENERATE_MAP_PLOTS)
        figure(2);
        subplot(3,1,1);
        plot(x_saved(k,1), x_saved(k,2), 'kx', x(1, :), x(2, :), 'r--', path(:,1), path(:,2), 'b:'); hold on; grid on;
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
end

f = figure();
set(f,'defaulttextInterpreter','latex') ;
set(f, 'defaultAxesTickLabelInterpreter','latex');  
set(f, 'defaultLegendInterpreter','latex');
set(f, "Position", [100,100,505,750]);
subplot(3,1,1);
% error plot
plot(Ts:Ts:N*Ts, error, 'DisplayName','Lane offset', 'LineWidth',1.5, 'Color','k');
grid on; hold on;
ylim([-1,1]); xlabel('time(s)'); ylabel('$\delta(t)$');
xline(0,'Alpha',0.8,'Color','k', 'HandleVisibility','off');
xlim([Ts,N*Ts]);
set(gca,'FontSize',18);
legend('Location','best', 'FontSize',18);
title("\textbf{Centerline driving}");

subplot(3,1,2);
% input plot
plot(Ts:Ts:N*Ts, u_saved, 'DisplayName', 'Steering angle', 'LineWidth',1.5, 'color', 'k');
grid on;
xlabel('time(s)'); ylabel('$\alpha_f(rad)$');
ylim([-0.05,0.05]); xlim([Ts,N*Ts]);
xline(0,'Alpha',0.8,'Color','k', 'HandleVisibility','off');
set(gca,'FontSize',18);
legend('Location','best', 'FontSize',18);

subplot(3,1,3);
% computational capacity plot
plot(Ts:Ts:N*Ts, compt_saved, 'DisplayName', 'Computational time', 'LineWidth',1, 'color', 'k');
grid on;
yline(mean(compt_saved), 'HandleVisibility','off', 'color', 'k', 'LineWidth', 2);
xlabel("time(s)"); ylabel("time steps(s)"); 
xlim([Ts,N*Ts]);
set(gca,'FontSize',18);
legend('Location','best', 'FontSize',18);

function f = objfunx2(U)
global x0 parameters pointsRefGlobal pointsOppositeLaneGlobal trafficDistance map priorPath path orientationRefVector yRefVector GENERATE_OPT_PLOTS partialCosts ay0 t0
    
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
    % output
    y = [0 1 0; 0 0 1]*xLocal';
    
    yref = [yRefVector'; orientationRefVector'];
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
    % input error
    % recalculation to lateral acceleration - ackermann angle
    Ay = parameters.vx^2*tan(U)/(parameters.Lr+parameters.Lf);
    f = f+(Ay/parameters.aymax*parameters.kappa(1)*Ay'/parameters.aymax)/size(y,2);
    
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