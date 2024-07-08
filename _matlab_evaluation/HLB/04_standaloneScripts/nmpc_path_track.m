clear; close all;

global modelID path
modelID = "eldm";

path = load("C:\git\hlb\Solutions\CruiseConcept\modelPredictiveControl\31_south_map.mat");
path = path.traj_local;
path = path(190:end,:);
orientation = diff(path(:,2))./diff(path(:,1));
orientation = [orientation; orientation(end)];
orientation = atan(movmean(orientation,25));
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

GENERATE_MAP_PLOTS = true;

possible_alfa_dy = 1;
possible_theta_dy = linspace(0,10,20);
possible_kappa = linspace(0,1,20);

global x0 parameters yRefVectorGlobal orientationRefVectorGlobal map priorPath 
parameters.Tsolve = 0.02; % only relevant when self-implemented solver is used
parameters.vx = 17; % m/s
parameters.Calfa_f = 10000; % * 180/pi(); % N/rad
parameters.Calfa_r = 10000; % * 180/pi(); % N/rad based on: https://www.google.com/url?sa=i&url=https%3A%2F%2Fwww.researchgate.net%2Ffigure%2FTire-Lateral-Force-vs-Slip-Angle-DOE-parameters_fig3_224323050&psig=AOvVaw2govbGsiWN6YjtVXFHzPqx&ust=1706344971038000&source=images&cd=vfe&opi=89978449&ved=0CBIQjRxqFwoTCNDd-8TU-oMDFQAAAAAdAAAAABAD

parameters.Lr = 1.8; % m
parameters.Lf = 1.2; % m
parameters.m = 1500; % kg
parameters.Iz = 2000; % kg*m^2

parameters.alfa = [2 0; 0 0];
parameters.kappa = 0.000;
parameters.beta = 0.0;
parameters.Sresidual = 0.0;

parameters.Th = 6;
parameters.Np = 20;
parameters.umax = 3*pi()/180; % absolute max road-wheel angle

parameters.costMode = "normal";
parameters.modelMode = "kinematic";

% planner
parameters.P_npDistances = [0.04, 0.116, 0.388];
parameters.P_ELDM = [0.39 0.24 0.12 -0.23 -0.39 -0.18 0.07; ...
    -0.39 -0.2 -0.1 0.43 0.6 0.18 0.07; ...
    0.08 0.04 0.07 -0.18 -0.18 0.03 0.07];
%parameters.P_ELDM = zeros(3,7);
    
Ts = 0.02;

U = zeros(1,parameters.Np);

T = 35;
N = T/Ts;
Td = parameters.Th/parameters.Np;
t_delaySteering = 2*Td;

options = optimoptions('fmincon','Display','iter','Algorithm','sqp');

M = 0;

if (parameters.modelMode == "dynamic")
    x_k = zeros(5,1);
    x_k(3) = path(1,1); x_k(4) = path(1,2); x_k(5) = path(1,3);
    vehicleState.X = x_k(3);
    vehicleState.Y = x_k(4);
    vehicleState.theta = x_k(5);
elseif (parameters.modelMode == "kinematic")
    x_k = zeros(3,1);
    x_k(1) = path(1,1); x_k(2) = path(1,2); x_k(3) = path(1,3);
    vehicleState.X = x_k(1);
    vehicleState.Y = x_k(2);
    vehicleState.theta = x_k(3);
end
x_k0 = x_k;
x_saved(1:2,:) = [x_k0'; x_k0'];
u_saved(1,1) = 0;
egoPath = [vehicleState.X vehicleState.Y];
globalPosteriorPath = egoPath;
% simulation loop
for k=1:N
    
    x0 = x_k;
    t_k = k*Ts;
    if (t_k/(M*Td) >= 1)
        % new plan - store previous global map points up until now
        if (k>1)
            nearestIndex = getNearestIndex(priorPath(:,1),priorPath(:,2), [vehicleState.X vehicleState.Y]);
            globalPosteriorPath = [globalPosteriorPath; priorPath(1:nearestIndex,:)];
        end
        % generate local map
        nearestIndex = getNearestIndex(path(:,1),path(:,2), [vehicleState.X vehicleState.Y]);
        scenario = functional_cutScenario(path, indexes, nearestIndex, vehicleState, [0,160]);
        priorPath = functional_planner (scenario, indexes, globalPosteriorPath, parameters);
        if (parameters.costMode == "costMap")
            map = generateMap(priorPath, [scenario(1:size(priorPath,1),indexes.corrX), scenario(1:size(priorPath,1),indexes.corrY)], [6, 2, -2, -2], 0.05, 0.1);
            map.Z(isnan(map.Z))=1;
        end
        %         
%         hold off;
%         h = pcolor(map.X,map.Y,map.Z);
%         clim([0 1]); 
%         set(h,'LineStyle','none');
%         view(2);
%         hold on;
%         plot(priorPath(:,1),priorPath(:,2), 'LineWidth',2,'color', 'r');
%         plot(scenario(:,indexes.corrX), scenario(:,indexes.corrY),'LineWidth',2,'LineStyle',':', 'color','k');
%         plot(vehicleState.X, vehicleState.Y, 'bo');
%         shg;
        % finding local path, interpolate to prediction horizon
        if (parameters.modelMode == "dynamic")
            distances = (sum((path(:,1:2)-[x_k(3) x_k(4)])'.^2)).^0.5;
        elseif (parameters.modelMode == "kinematic")
            distances = (sum((path(:,1:2)-[x_k(1) x_k(2)])'.^2)).^0.5;
        end
        nearestIndex = find(distances==min(distances),1);
        if (~isempty(nearestIndex))
            % there is a valid point found globally, transform trajectory to
            % local
            dx = parameters.vx*parameters.Th/parameters.Np;
            xend = dx*parameters.Np; % distance during one optimization cycle
            if (parameters.modelMode == "dynamic")
                T = [cos(x_k(5)) sin(x_k(5)); -sin(x_k(5)) cos(x_k(5))];
                P = [x_k(3) x_k(4)];
            else
                T = [cos(x_k(3)) sin(x_k(3)); -sin(x_k(3)) cos(x_k(3))];
                P = [x_k(1) x_k(2)];
            end
            localPath = (priorPath(:,1:2)-P)*T';
            localDistances = (sum(localPath'.^2)).^0.5;
            endPoint = find(localDistances(1:end)>=xend,1);
            xRefVector = linspace(dx,xend, parameters.Np);
            yRefVector = spline(localPath(1:endPoint,1), localPath(1:endPoint,2), xRefVector);
            orientationRefVector = movmean(diff(yRefVector)./diff(xRefVector), 10);
            orientationRefVector = [orientationRefVector orientationRefVector(end)];
            % transform back to global
            pointsRefGlobal = [xRefVector' yRefVector']*T+P;
            yRefVectorGlobal = pointsRefGlobal(:,2);
            if (parameters.modelMode == "dynamic")
                orientationRefVectorGlobal = orientationRefVector'+x_k(5);
            else
                orientationRefVectorGlobal = orientationRefVector'+x_k(3);
            end
        end
        U = zeros(1,parameters.Np);
        M = M+1;
        %[U, fval(M)] = fmincon(@objfunx, U, [], [], [], [], -parameters.umax*ones(1,parameters.Np), ones(1,parameters.Np)*parameters.umax);
        for o=1:1
            [U, fval(M), exitFlag] = fminsearch(@objfunx2, U); 
        end
        x_k0 = x_k;
        for i=1:parameters.Np
            t_i = (i-1)*Td;
            u_i = mu(U, t_i);
            if (i==1)
                x(:,i) = kinematicSingleTrack(u_i, x_k, Td);
            else
                x(:,i) = kinematicSingleTrack(u_i, x(:,i-1), Td);
            end
            %x(:, i) = phi(t_i, U, x0);
        end 
    end
    u_k = mu(U, t_k-(M-1)*Td); % receeding horizon, shifting with dela
    %u_k = 0.5*(mu(U, t_k-(M-1)*Td+Td)+mu(U, t_k-(M-1)*Td+2*Td));
    % saving for debugging
    u_saved(k) = u_k;
    x_saved(k,1:size(x_k,1)) = x_k';
    th_saved(k) = t_k-(M-1)*Td;

    % state update
    %x_k = F(x_k, u_k, Ts); 
    x_k = kinematicSingleTrack(u_k, x_k, Ts);
    if (parameters.modelMode == "dynamic")
        T = [cos(x_k(5)) sin(x_k(5)); -sin(x_k(5)) cos(x_k(5))];
    else
        T = [cos(x_k(3)) sin(x_k(3)); -sin(x_k(3)) cos(x_k(3))];
    end
    localPriorPath = (priorPath-[vehicleState.X vehicleState.Y])*T';
    error(k) = spline(localPriorPath(:,1), localPriorPath(:,2), 0);
    
    if(GENERATE_MAP_PLOTS)
        subplot(3,1,1);
%         hold off;
%         h = pcolor(map.X,map.Y,map.Z);
%         clim([0 1]); 
%         set(h,'LineStyle','none');
%         view(2);
%         hold on;
        if (parameters.modelMode == "dynamic")
            plot(x_saved(k,3), x_saved(k,4), 'kx', x(3, :), x(4, :), 'r--', priorPath(:,1), priorPath(:,2), 'b:', globalPosteriorPath(:,1), globalPosteriorPath(:,2), 'ko'); hold on; grid on;
        else
            plot(x_saved(k,1), x_saved(k,2), 'kx', x(1, :), x(2, :), 'r--', priorPath(:,1), priorPath(:,2), 'b:', globalPosteriorPath(:,1), globalPosteriorPath(:,2), 'ko'); hold on; grid on;
        end
            %xlim([min(x_saved(:,3)) max(10,max(x(:,3)))]);
        title(strcat("error is:", num2str(error(k))));
        subplot(3,1,2);
        if (parameters.modelMode == "dynamic")
            plot(x_saved(k,3), u_saved(k), 'kx', x(3, :), U, 'r--'); hold on; grid on; xlabel('X(m)'); ylabel('\delta_f')
        else
            plot(x_saved(k,1), u_saved(k), 'kx', x(3, :), U, 'r--'); hold on; grid on; xlabel('X(m)'); ylabel('\delta_f')
        end
        subplot(3,1,3);
        if (parameters.modelMode == "dynamic")
            plot(x_saved(k,3), error(k), 'kx'); grid on; hold on; ylabel('error(m)');
        else
            plot(x_saved(k,1), error(k), 'kx'); grid on; hold on; ylabel('error(m)');
        end
        pause(0.1);
    end

    % update vehicle state struct
    if (parameters.modelMode == "dynamic")
        vehicleState.X = x_k(3);
        vehicleState.Y = x_k(4);
        vehicleState.theta = x_k(5);
    elseif (parameters.modelMode == "kinematic")
        vehicleState.X = x_k(1);
        vehicleState.Y = x_k(2);
        vehicleState.theta = x_k(3);
    end
    egoPath = [egoPath; [vehicleState.X vehicleState.Y]];

    clc;
    fprintf("Iteration number %d / %d\n", k, N);
end

function f = objfunx(U)
    global x0 parameters yRefVectorGlobal map priorPath path
    f = 0;
    t_k1 = 0;
    Td = parameters.Th/parameters.Np;

    for k=1:parameters.Np
        t_k = k*Td;
        u_k = mu(U, t_k);
        f = f+l(phi(t_k, U, x0), u_k, t_k)*(t_k-t_k1);
        t_k1 = t_k;
        x(:, k) = phi(t_k, U, x0);
    end 
    f = f/parameters.Th;
    f = f + S(phi(parameters.Th,U,x0), parameters.Th);
    
%     if(parameters.costMode=="costMap")
%         hold off;
%         h = pcolor(map.X,map.Y,map.Z);
%         clim([0 1]); 
%         set(h,'LineStyle','none');
%         view(2);
%         hold on;
%     end
%     if (parameters.modelMode == "dynamic")
%         plot(x(3,:), x(4,:), 'r', x(3,:), yRefVectorGlobal(1:end), 'k', priorPath(:,1), priorPath(:,2), 'g', path(:,1), path(:,2), 'k:', 'LineWidth',2);
%     else
%         plot(x(1,:), x(2,:), 'r', x(1,:), yRefVectorGlobal(1:end), 'k');
%     end
%     title(num2str(f));
%     xlim([min(priorPath(:,1))-3, max(priorPath(:,1))+3]);
%     ylim([min(priorPath(:,2))-3, max(priorPath(:,2))+3]);
%     pause(0.05);
end

function f = objfunx2(U)
global x0 parameters yRefVectorGlobal map priorPath path orientationRefVectorGlobal
    
    Td = parameters.Th/parameters.Np;
    x_k = x0;
    for k=1:parameters.Np
        t_k = k*Td;
        u_k = mu(U, t_k);
        x_k = kinematicSingleTrack(u_k, x_k, Td);
        x(:, k) = x_k;
    end 
    % output
    y = [0 1 0; 0 0 1]*x;
    % cost 
    yref = [yRefVectorGlobal'; orientationRefVectorGlobal'];
    f = 0;
    for n=1:size(y,1)
        yn = y(n,:); yrefn = yref(n,:);
        f = f+(yn-yrefn)*parameters.alfa(n,n)*(yn-yrefn)';
    end
    f = f+U*parameters.kappa*U';
    
    if(parameters.costMode=="costMap")
        hold off;
        h = pcolor(map.X,map.Y,map.Z);
        clim([0 1]); 
        set(h,'LineStyle','none');
        view(2);
        hold on;
    end
    if (parameters.modelMode == "dynamic")
        plot(x(3,:), x(4,:), 'r', x(3,:), yRefVectorGlobal(1:end), 'k', priorPath(:,1), priorPath(:,2), 'g', path(:,1), path(:,2), 'k:', 'LineWidth',2);
    else
        plot(x(1,:), x(2,:), 'r', x(1,:), yRefVectorGlobal(1:end), 'k', priorPath(:,1), priorPath(:,2), 'g', path(:,1), path(:,2), 'k:', 'LineWidth',2);
    end
    title(num2str(f));
    xlim([min(x(1,:))-3, max(x(1,:))+3]);
    ylim([min(x(2,:))-1, max(x(2,:))+1]);
    pause(0.05);

end

function fl = l(x,u, t_k)
    global parameters map
    y = output(x);
    if (parameters.costMode == "normal")
        y_ref = calculateReference(t_k);
        fl = (y-y_ref)'*parameters.alfa*(y-y_ref)+parameters.kappa*u^2;
    elseif (parameters.costMode == "costMap")
        % get nearest grid point
        if (parameters.modelMode == "dynamic")
            distanceX = ((x(3)-map.X(1,:)).^2).^0.5;
            idxX = find(distanceX==min(distanceX),1);
            distanceY = ((x(4)-map.Y(:,idxX)).^2).^0.5;        
            idxY = find(distanceY==min(distanceY),1);
        else
            distanceX = ((x(1)-map.X(1,:)).^2).^0.5;
            idxX = find(distanceX==min(distanceX),1);
            distanceY = ((x(2)-map.Y(:,idxX)).^2).^0.5;        
            idxY = find(distanceY==min(distanceY),1);
        end
        fl = parameters.alfa(1,1)*map.Z(idxY,idxX)+parameters.kappa*(u./parameters.umax)^2;
    else
        fl = 0;
    end
end

function y = output(x)
global parameters
    if (parameters.modelMode == "dynamic")
        y(1,1) = x(4,1);
        y(2,1) = x(5,1);
    elseif (parameters.modelMode == "kinematic")
        y(1,1) = x(2,1);
        y(2,1) = x(3,1);
    end
end

function s = S(x,T)
    global parameters
    y = output(x);
    y_ref = calculateReference(T);
    s = parameters.Sresidual + parameters.beta*(y-y_ref)'*(y-y_ref);
end

function phi_tk = phi(t_N, U, x0)
    % ODE solution function
    global parameters
    N = t_N/parameters.Tsolve;
    x_k = x0; 
    if (N==0)
        phi_tk = x_k;
    else
%         for k=0:(N-1)
%             u_k = mu(U, k*parameters.Tsolve);
%             t_k = (k+1)*parameters.Tsolve;
%             x_k = F(x_k, u_k, parameters.Tsolve);
%         end
        phi = ode45(@(t_N, x_k)ODE(t_N, x_k, U), [0, t_N], x0);
        phi_tk = phi.y(:,end);
        %phi_tk = x_k;
    end
end

function y_err = calculateReference(t_k)
    global parameters yRefVectorGlobal orientationRefVectorGlobal
    % reference point calculation
    tPredStep = parameters.Th/parameters.Np;
    n_act = floor(t_k/tPredStep);
    if (n_act >= parameters.Np)
        y_refCalculated = yRefVectorGlobal(end);
        orientation_refCalculated = orientationRefVectorGlobal(end);
    else
    % n_act+1 in indexing is needed, as time starts from 0
    y_refCalculated = interp1([n_act*tPredStep, (n_act+1)*tPredStep], ...
        [yRefVectorGlobal(n_act),yRefVectorGlobal(n_act+1)], ...
        t_k);
    orientation_refCalculated = interp1([n_act*tPredStep, (n_act+1)*tPredStep], ...
        [orientationRefVectorGlobal(n_act),orientationRefVectorGlobal(n_act+1)], ...
        t_k);
    end
    y_err = [y_refCalculated; orientation_refCalculated];
end

function u_k = mu(U, t_k)
    global parameters
    Td = parameters.Th/parameters.Np;
    N = max(1,floor(t_k/Td));
%     if (N<size(U,2))
%         % linear interpolation is possible
%         if (t_k < N*Td)
%             u_k = (U(N+1)-U(N))/Td*t_k+U(N);
%         else
%             u_k = (U(N+1)-U(N))/Td*(t_k-N*Td)+U(N);
%         end
%     else
        u_k = U(N);
    %end
end

function x_k = F(x_k1, U, Td)
    global parameters
    dxdt = ODE(Td, x_k1, U);
    x_k = x_k1+dxdt*Td;
end

function dxdt = ODE(t, x, U)
    global parameters
    u = mu(U, t);
    if (parameters.modelMode == "dynamic")
        A = - (parameters.Calfa_f*cos(u(1))+parameters.Calfa_r)/(parameters.m*parameters.vx);
        B = (-parameters.Lf*parameters.Calfa_f*cos(u(1))+parameters.Lr*parameters.Calfa_r)/(parameters.m*parameters.vx) - parameters.vx;
        C = (-parameters.Lf*parameters.Calfa_f*cos(u(1))+parameters.Lr*parameters.Calfa_r)/(parameters.Iz*parameters.vx);
        D = (-parameters.Lf^2*parameters.Calfa_f*cos(u(1))+parameters.Lr^2*parameters.Calfa_r)/(parameters.Iz*parameters.vx);
        E = (parameters.Calfa_f*cos(u(1)))/parameters.m;
        F = (parameters.Lf*parameters.Calfa_f*cos(u(1)))/parameters.Iz;
    
        dxdt(1,1) = A*x(1)+C*x(2)+E*u;
        dxdt(2,1) = B*x(1)+D*x(2)+F*u;
        dxdt(3,1) = parameters.vx*cos(x(5))-x(1)*sin(x(5));
        dxdt(4,1) = parameters.vx*sin(x(5))+x(1)*cos(x(5));
        dxdt(5,1) = x(2);
    elseif (parameters.modelMode == "kinematic")
        % x = [X,Y,theta]
        dxdt(1,1) = cos(x(3))*parameters.vx;
        dxdt(2,1) = sin(x(3))*parameters.vx;
        dxdt(3,1) = parameters.vx*tan(u(1))/(parameters.Lf+parameters.Lr);
    end
end

function [c, ceq] = g(U)

% ODE solution function
global x0 parameters

Td = parameters.Th/parameters.Np;

for k=1:parameters.Np
    t_k = k*Td;
    x_k = phi(t_k, U, x0);
    c(k) = x_k(4,1)-1.1;
end

ceq = [];
end

function map = generateMap(priorPath, centerLine, boundaries, gridStepMain, kernelStep)
% cost generation: from map with gaussian distribution: f(x) =
% a*exp(-(x-b)^2/2c^2)
% where: a= 1/(sigma*sqrt(2*pi)) b=mu, c=sigma
routeLength = sum((diff(priorPath(:,1)).^2+diff(priorPath(:,2)).^2).^0.5);
N = floor(routeLength/gridStepMain);
local_theta = atan(diff(priorPath(1:2,2))/diff(priorPath(1:2,1)));
resampledPath(1,:) = priorPath(1,:);
for i=2:N
    dx = gridStepMain*cos(local_theta);
    if (i==2)
        resampledPath(i,1) = priorPath(1,1)+dx;
        resampledPath(i,2) = spline(priorPath(1:2,1), priorPath(1:2,2), resampledPath(i,1));
    else
        resampledPath(i,1) = resampledPath(i-1,1)+dx;
        resampledPath(i,2) = spline(priorPath(:,1), priorPath(:,2), resampledPath(i,1));
    end
    local_theta = atan(diff(resampledPath(i-1:i,2))/diff(resampledPath(i-1:i,1)));
end
resampledCenterLine(:,1) = resampledPath(:,1);
resampledCenterLine(:,2) = spline(centerLine(:,1), centerLine(:,2), resampledPath(:,1));

% calculate boundary based cost in the Frenét frame of the centerline first
% kernel function shape: inverse gaussian shape
wDrivable = 0.0; wEgo = 0.0; wInstinctive = 1.0;
for i=1:N
   
    x = boundaries(4):kernelStep:boundaries(1);
    % drivable surface
    x_normalized = x/max(abs(x));
    costDrivableSurface = kernelInverseGaussian(x_normalized, 0.25, boundaries(2)/max(abs(x)));
    % ego lane - without leaving it
    costEgoLane = kernelInverseGaussian(x_normalized, 0.125, 0);
    % instinctive path distance cost
    costInstinctivePath = kernelInverseGaussian(x_normalized, 0.5, (resampledPath(i,2)-resampledCenterLine(i,2))/max(abs(x)));
    % weighted sum
    cost = wDrivable*costDrivableSurface+wEgo*costEgoLane+wInstinctive*costInstinctivePath;
    if(i==1)
        local_theta = atan(diff(resampledCenterLine(1:2,2))/diff(resampledCenterLine(1:2,1)));
        mapV = [resampledCenterLine(i,1)-sin(local_theta)*x' resampledCenterLine(i,2)+cos(local_theta)*x', cost'];
    else
        local_theta = atan(diff(resampledCenterLine(i-1:i,2))/diff(resampledCenterLine(i-1:i,1)));
        mapV = [mapV; ...
            [resampledCenterLine(i,1)-sin(local_theta)*x' resampledCenterLine(i,2)+cos(local_theta)*x', cost']];
    end
end
% interpolating grid 
xv = linspace(min(mapV(:,1)), max(mapV(:,1)), N);
yv = linspace(min(mapV(:,2)), max(mapV(:,2)), N);
[X,Y] = meshgrid(xv, yv);
Z = griddata(mapV(:,1),mapV(:,2),mapV(:,3),X,Y);
map.X = X; map.Y = Y; map.Z=Z;
end

function cost = kernelInverseGaussian(x, sigma, mu)
    a = 1/(sigma*sqrt(2*pi()));
    b = mu;
    c = sigma;
    cost = 1-a*exp(-(x-b).^2/(2*c^2));
    cost = cost/(1-min(cost));
    cost = cost-min(cost); % shift to zero local extremum    
end

function cost = kernelGaussian(x, sigma, mu)
    a = 1/(sigma*sqrt(2*pi()));
    b = mu;
    c = sigma;
    cost = a*exp(-(x-b).^2/(2*c^2));
end

function U = lmpc(pathLocal, pathOrientation, Np, Nc, rw, q, L, Ts,  x_k1, xa, delta0, u_k1)
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
    Ad = [1 0 Ts 0 0; 0 1 0 Ts 0; 0 0 1 0 -Ts/L*x_k1(3)^2*tan(delta0); 0 0 0 1 0; 0 0 0 0 1];
    Bd = [0 0 0 Ts/L*x_k1(3)^2*1/(cos(delta0))^2 Ts/L*x_k1(3)*1/(cos(delta0))^2]';
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
    U(1) = u_k1+dU(1);
    for i=2:length(dU)
        U(i) = U(i-1)+dU(i);
    end
end
end

function x_k = kinematicSingleTrack(u_k, x_k1, Ts)
global parameters
x_k(1,1) = x_k1(1,1)+Ts*parameters.vx*cos(x_k1(3));
x_k(2,1) = x_k1(2,1)+Ts*parameters.vx*sin(x_k1(3));
x_k(3,1) = x_k1(3,1)+parameters.vx*Ts*tan(u_k(1))/(parameters.Lr+parameters.Lf);
end
