clear; close all;
global x0 parameters y_ref
parameters.Tsolve = 0.02; % only relevant when self-implemented solver is used
parameters.vx = 30; % m/s
parameters.Calfa_f = 100000; % * 180/pi(); % N/rad
parameters.Calfa_r = 100000; % * 180/pi(); % N/rad based on: https://www.google.com/url?sa=i&url=https%3A%2F%2Fwww.researchgate.net%2Ffigure%2FTire-Lateral-Force-vs-Slip-Angle-DOE-parameters_fig3_224323050&psig=AOvVaw2govbGsiWN6YjtVXFHzPqx&ust=1706344971038000&source=images&cd=vfe&opi=89978449&ved=0CBIQjRxqFwoTCNDd-8TU-oMDFQAAAAAdAAAAABAD

parameters.Lr = 1.8; % m
parameters.Lf = 1.2; % m
parameters.m = 1500; % kg
parameters.Iz = 2000; % kg*m^2

parameters.alfa = [0.1 0; 0 0.001];
parameters.kappa = 1;
parameters.beta = 0.5;

parameters.Th = 6;
parameters.Np = 10;

y_ref = [1; 0];
    
Ts = 0.02;
U = zeros(1,parameters.Np);

T = 20;
N = T/Ts;
Td = parameters.Th/parameters.Np;
t_delaySteering = 0;

options = optimoptions('fmincon','Display','iter','Algorithm','sqp');

x_k = zeros(5,1);
x_k0 = x_k;
M = 0;

for k=1:N
    x0 = x_k;
    t_k = k*Ts;
    if (t_k/(M*Td) >= 1)
        U = zeros(1,parameters.Np);
        M = M+1;
        [U, fval(M)] = fmincon(@objfunx, U, [], [], [], [], -0.3*ones(1,parameters.Np), ones(1,parameters.Np)*0.3);
        x_k0 = x_k;
        for i=1:parameters.Np
            t_i = (i-1)*Td;
            u_i = mu(U, t_i);
            x(:, i) = phi(t_i, U, x0);
        end 
    end
    u_k = mu(U, t_k-(M-1)*Td+t_delaySteering); % receeding horizon, shifting with delay
    % saving for debugging
    u_saved(k) = u_k;
    x_saved(k,1:5) = x_k';
    th_saved(k) = t_k-(M-1)*Td;

    % state update
    x_k = F(x_k, u_k, Ts);  
    subplot(2,1,1);
    plot(x_saved(k,3), x_saved(k,4), 'kx', x(3, :), x(4, :), 'r--'); hold on; grid on;
    subplot(2,1,2);
    plot(x_saved(k,3), u_saved(k), 'kx', x(3, :), U, 'r--'); hold on; grid on; xlabel('X(m)'); ylabel('\delta_f')
    pause(0.1);

    clc;
    fprintf("Iteration number %d / %d", k, N);
end



function f = objfunx(U)
    global x0 parameters
    f = 0;
    t_k1 = 0;
    Td = parameters.Th/parameters.Np;

    for k=1:parameters.Np
        t_k = (k-1)*Td;
        u_k = mu(U, t_k);
        f = f+l(phi(t_k, U, x0), u_k, t_k)*(t_k-t_k1);
        t_k1 = t_k;
        x(:, k) = phi(t_k, U, x0);
    end    
    f = f + S(phi(parameters.Th,U,x0), parameters.Th);
end

function fl = l(x,u, t_k)
    global y_ref parameters
    y = output(x);
    fl = (y-y_ref)'*parameters.alfa*(y-y_ref)+parameters.kappa*u^2;
end

function y = output(x)
    y(1,1) = x(4,1);
    y(2,1) = x(5,1);
end

function s = S(x,T)
    global parameters y_ref
    y = output(x);
    s = parameters.beta*(y-y_ref)'*(y-y_ref);
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

function u_k = mu(U, t_k)
    global parameters
    Td = parameters.Th/parameters.Np;
    N = max(1,floor(t_k/Td));
    if (N<size(U,2))
        % linear interpolation is possible
        if (t_k < N*Td)
            u_k = (U(N+1)-U(N))/Td*t_k+U(N);
        else
            u_k = (U(N+1)-U(N))/Td*(t_k-N*Td)+U(N);
        end
    else
        u_k = U(N);
    end
end

function x_k = F(x_k1, U, Td)
    global parameters
    dxdt = ODE(Td, x_k1, U);
    x_k = x_k1+dxdt*Td;
end

function dxdt = ODE(t, x, U)
    global parameters
    u = mu(U, t);
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