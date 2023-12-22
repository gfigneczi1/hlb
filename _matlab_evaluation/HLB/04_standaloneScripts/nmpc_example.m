clear; close all;
global x0 parameters x_ref
parameters.Tsolve = 0.02;
parameters.v = 0.1;
parameters.beta = 0;
parameters.R = 0.8;
parameters.alfa = 1;
parameters.kappa = 0.08;
parameters.Th = 12;
parameters.Np = 15;
parameters.B = 1;
parameters.H = 0.18;
parameters.sigma_c0 = 0.3;
parameters.W = 0.25;
parameters.gamma = 0.5;

x_ref = [0.4; 0.6];
    
Ts = 0.1;
U = zeros(1,15);

T = 40;
N = T/Ts;
Td = 0.8;
Tdelay = 20;
D = Tdelay/Ts;

options = optimoptions('fmincon','Display','iter','Algorithm','sqp');

x_k = [0.7; 0.5];
x_k0 = x_k;
M = 0;

for k=1:(N+D)
    x0 = x_k;
    t_k = k*Ts;
    if (t_k/(M*Td) >= 1 && k>=D)
        M = M+1;
        %[U, fval(M)] = fminsearch(@objfunx, U);
        [U, fval(M)] = fmincon(@objfunx, U, [], [], [], [], zeros(1,15), ones(1,15)*0.3, @g);
        x_k0 = x_k;
    elseif(t_k/(M*Td) >= 1 && k<D)
        M = M+1;
        U = zeros(1,15);
        x_k0 = x_k;
    end 
    u_k = mu(U, t_k-(M-1)*Td); % receeding horizon
    % saving for debugging
    u_saved(k) = u_k;
    x_saved(k,1:2) = x_k';
    th_saved(k) = t_k-(M-1)*Td;

    % state update
    x_k = F(x_k, u_k, Ts);  

    clc;
    fprintf("Iteration number %d / %d", k, N+D);
end



function f = objfunx(U)
    global x0 parameters
    f = 0;
    t_k1 = 0;
    Td = parameters.Th/parameters.Np;

    for k=1:parameters.Np
        t_k = k*Td;
        u_k = mu(U, t_k);
        f = f+l(phi(t_k, U, x0), u_k, t_k)*(t_k-t_k1);
        t_k1 = t_k;
    end
    
    f = f + S(phi(parameters.Th,U,x0), parameters.Th);
end

function fl = l(x,u, t_k)
    global x_ref parameters
    fl = parameters.alfa*(x-x_ref)'*(x-x_ref)+parameters.kappa*u^2;
end

function s = S(x,T)
    global parameters x_ref
    s = parameters.R*parameters.v^2+parameters.beta*(x-x_ref)'*(x-x_ref);
end

function phi_tk = phi(t_N, U, x0)
    % ODE solution function
    global parameters
    N = t_N/parameters.Tsolve;
    x_k = x0; 
    if (N==0)
        phi_tk = x_k;
    else
        for k=0:(N-1)
            u_k = mu(U, k*parameters.Tsolve);
            t_k = (k+1)*parameters.Tsolve;
            x_k = F(x_k, u_k, parameters.Tsolve);
        end
        %phi = ode45(@(t_N, x_k)ODE(t_N, x_k, U), [0, t_N], x0);
        %phi_tk = phi.y(:,end);
        phi_tk = x_k;
    end
end

function u_k = mu(U, t_k)
    global parameters
    Td = parameters.Th/parameters.Np;
    N = max(1,floor(t_k/Td));
    u_k = U(N);
end

function x_k = F(x_k1, U, Td)
    global parameters
    dxdt = ODE(Td, x_k1, U);
    x_k = x_k1+dxdt*Td;
end

function xe = sigma_e(x)
    global parameters
    xe = parameters.sigma_c0+parameters.H*(1+1.5*(x/parameters.W-1)-0.5*(x/parameters.W-1)^3);
end

function x_out = fi(x)
    global parameters
    x_out = parameters.gamma*sign(x)*sqrt(abs(x));
end

function dxdt = ODE(t, x, U)
    global parameters
    u = mu(U, t);
    dxdt(1,1) = parameters.B*(sigma_e(x(1))-x(2)-u);
    dxdt(2,1) = 1/parameters.B*(x(1)-fi(x(2)));
end

function [c, ceq] = g(U)

% ODE solution function
global x0 parameters

Td = parameters.Th/parameters.Np;

for k=1:parameters.Np
    t_k = k*Td;
    x_k = phi(t_k, U, x0);
    c(k) = 0.4-parameters.v-x_k(2);
end

ceq = [];
end