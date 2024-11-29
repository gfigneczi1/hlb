function [u,f_value] = mpcCtrl_trans_GP(x,xrot,t,para,i)

global mpcobj 
mpcobj.A = para.trans.A;
mpcobj.B = para.trans.B;
mpcobj.Q = para.trans.Q;
mpcobj.R = para.trans.R;
mpcobj.P = para.trans.P;
mpcobj.H = para.trans.H;
% mpcobj.H = H;
mpcobj.dT =  para.dT;
mpcobj.i = i;
mpcobj.N = para.N;

mpcobj.dynmodel = para.gpmodel_trans;
mpcobj.xrot = xrot;


mpcobj.path = para.path;
mpcobj.n = 6;  % state dimension
mpcobj.m = 3;  % input dimension 


mpcobj.lb = -inf;  % input constraint, low 
mpcobj.ub = inf;   % input constraint, up

[u,f_value] = mpcSolver_trans_GP(x,t);


end