
% global m g Ct Cq l Cw e 


%% simulation params
para.T  = 50;       % simulation time (seconds)
para.dT = 0.05;      % step size (seconds)
para.time = 0:para.dT:para.T;  % time span
para.N = length(para.time);
% 
% para.rot.sampleInt   = para.dT;     % inner loop faster
% para.trans.sampleInt = 10*para.dT;  % outer loop

%% quad dyanmics params
para.m  = 1.9;  % kg, the mass of quad
para.g  = 9.81; % kg m^2, the gravity
para.Jx = 0.0059;
para.Jy = 0.0059;
para.Jz = 0.0107; % inertia moments
% para.J = [Jx Jy Jz];

% para.L  = 0.315; % rotor-to-rotor distance
% para.l  = para.L/sqrt(2);

para.l  = 0.25;  % m, 90 deg cross configuration

para.e  = 0;   % >0, rotors are above the center of mass;
          % <0, rotors are below the center of mass.
          
para.Ct = 1e-05; % constant cofficient about force & rotor speed w
para.Cq = 1e-06; % constant cofficient about torque & rotor speed w  
para.Cw = 0.1924; %rho * cd * SA * (1/2)

%% MPC settings
% para.rot.H   = 10;
% para.trans.H = 10;

para.rot.H   = 1;
para.trans.H = 5;

para.trans.Q = diag([1 1 20 1 1 20]);
para.trans.R = diag([1 1 1]);

para.rot.Q = diag([1 1 1 1 1 1]);
para.rot.R = diag([1 1 1]);

para.opts = optimset('Display','off');

