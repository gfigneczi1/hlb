function [u,f_value] = mpcCtrl_trans_lGP(x,xrot,t,para,i)

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

mpcobj.model_x = para.model_x;
mpcobj.model_y = para.model_y;
mpcobj.model_z = para.model_z;

% mpcobj.Smodel_x = para.Smodel_x;
% mpcobj.Smodel_y = para.Smodel_y;
% mpcobj.Smodel_z = para.Smodel_z;
% 

mpcobj.xrot = xrot;

mpcobj.path = para.path;
mpcobj.n = 6;  % state dimension
mpcobj.m = 3;  % input dimension 


mpcobj.lb = -inf;  % input constraint, low 
mpcobj.ub = inf;   % input constraint, up

mpcobj.con = [-3 -3 0 -2 -2 -2 3 3 4 2 2 2];  % min (x y z vx vy vz),max(x y z vx vy vz)


[u,f_value] = mpcSolver_trans_lGP(x,t);


end