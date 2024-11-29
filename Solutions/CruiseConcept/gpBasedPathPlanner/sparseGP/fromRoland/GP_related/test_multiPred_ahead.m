clc;clear; close all

load gpmodel_init.mat

set_parameters;
init_LTI;

xrot = rand(6,1);
%% Dyn settings

para.trans.A  = LTI_trans.A;
para.trans.B  = LTI_trans.B;  % discrete 

para.rot.A  =  LTI_rot.A;
para.rot.B  =  LTI_rot.B;  % discrete 

para.model_x = model_x;
para.model_y = model_y;
para.model_z = model_z;
para.Smodel_x = Smodel_x;
para.Smodel_y = Smodel_y;
para.Smodel_z = Smodel_z;


%%
mpcobj.A = para.trans.A;
mpcobj.B = para.trans.B;
mpcobj.Q = para.trans.Q;
mpcobj.R = para.trans.R;
mpcobj.P = para.trans.P;
mpcobj.H = para.trans.H;
mpcobj.H = 10;
mpcobj.dT =  para.dT;
% mpcobj.i = i;
mpcobj.N = para.N;
mpcobj.xrot = xrot;

mpcobj.model_x = para.model_x;
mpcobj.model_y = para.model_y;
mpcobj.model_z = para.model_z;

mpcobj.Smodel_x = para.Smodel_x;
mpcobj.Smodel_y = para.Smodel_y;
mpcobj.Smodel_z = para.Smodel_z;

mpcobj.n = 6;  % state dimension
mpcobj.m = 3;  % input dimension 

%%
x0 = rand(6,1);
mu_x0 = x0;
var_x0 = zeros(mpcobj.n);
uk = rand(3,1);
uk = repmat(uk,1,mpcobj.H);
t  = 0;

[mu_xk,var_xk] = multi_pred_ahead(mu_x0,var_x0,uk,t,mpcobj);

var_x = squeeze(var_xk(1,1,:));
std_x = sqrt(var_x);

time  =  1:mpcobj.H+1;

fillColor = [236 241 249]/255;

% 1-D
figure;
hold on
fill([time'; time(end:-1:1)'], [mu_xk(1,:)';mu_xk(1,end:-1:1)']...
               + 2*[std_x; -std_x(end:-1:1)], fillColor,'EdgeColor','none');
plot(time,mu_xk(1,:),'linewidth',2); 

% 2-D
phi = 0:pi/20:2*pi;
R = 0.1;
figure;
hold on
for i = 1:mpcobj.H+1
     R = 2*std_x(i);
     o_x = R*cos(phi)+mu_xk(1,i);
     o_y = R*sin(phi)+mu_xk(2,i);
     
     fill(o_x,o_y,fillColor,'EdgeColor','none');
end
plot(mu_xk(1,:),mu_xk(2,:),'linewidth',2); 
plot(mu_xk(1,1),mu_xk(2,1),'ro','linewidth',2,'markersize',10); 

