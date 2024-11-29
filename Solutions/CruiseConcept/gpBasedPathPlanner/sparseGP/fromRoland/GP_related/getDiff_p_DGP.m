function training_out_p = getDiff_p_DGP(x,u,mu_sp_hist,state_dot_hist,para,LTI)

% x,y,z, vx,   vy,   vz, phi, theta, psi,  p,  q,  r
% 1 2 3   4     5     6   7     8    9     10  11  12

% state_hist-----  N*12 
% statedot_hist-----  N*12
% F_hist --------  1*N

N  = size(state_dot_hist,2);

observations = state_dot_hist(4:6,:)' ;   % acceleration in x,y,z direction Nx3

% phi = state_hist(:,4);
% theta = state_hist(:,5);
% psi  = state_hist(:,6);
% 
% phi  = res.x.rot(1,1:end-1)';
% theta =  res.x.rot(2,1:end-1)';
% F_hist = res.u.trans(3,:)';
% 
% nominal_x = -para.g*theta;
% nominal_y = para.g*phi;
% nominal_z = -F_hist/para.m+para.g;
H = para.trans.H;

xtrans = x(1:6,1:end-H);
utrans = u(1:3,:);

dx_trans = LTI.trans.A*xtrans + LTI.trans.B*utrans;
nominal = dx_trans(4:6,:)';

% 
% nominal_x = -F_hist'.*R13/m;
% nominal_y = -F_hist'.*R23/m;
% nominal_z = -F_hist'.*R33/m + g;

diff_x = observations(:,1) - nominal(:,1);
diff_y = observations(:,2) - nominal(:,2);
diff_z = observations(:,3) - nominal(:,3);

training_out_p = [diff_x diff_y diff_z];   % N*3
