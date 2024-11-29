function [u,f_value] = mpcCtrl_trans_noComp(x,t,para,i)

% ref_sampled = interp1(para.time,xd',t+(0:para.H)*para.dT)'; % 2x(H+1)

global mpcobj
mpcobj.A = para.trans.A;
mpcobj.B = para.trans.B;
mpcobj.Q = para.trans.Q;
mpcobj.R = para.trans.R;
mpcobj.P = para.trans.P;
mpcobj.H = para.trans.H;
mpcobj.dT =  para.dT;
mpcobj.i = i;
mpcobj.N = para.N;

mpcobj.path = para.path;
mpcobj.n = 6;  % state dimension
mpcobj.m = 3;  % input dimension 

mpcobj.lb = -inf;  % input constraint, low 
mpcobj.ub = inf;   % input constraint, up

[u,f_value] = mpcSolver_trans_noComp(x,t);


end