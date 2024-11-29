function training_out_a = getDiffOnline_a(x,u,state_dot_hist,para)

% x,y,z, vx,   vy,   vz, phi, theta, psi,  p,  q,  r
% 1 2 3   4     5     6   7     8    9     10  11  12


% state_hist-----  N*12
% statedot_hist-----  N*12


Jx = para.Jx;
Jy = para.Jy;
Jz = para.Jz;

tau_x = u(4,:);
tau_y = u(5,:);
tau_z = u(6,:);

p    = x(10,:);
q    = x(11,:);
r    = x(12,:);

observations = state_dot_hist(10:12,:)' + 1e-8*randn(1,3);   % 

Jpdot = tau_x - (Jz-Jy) * q.*r;
Jqdot = tau_y - (Jx-Jz) * p.*r;
Jrdot = tau_z - (Jy-Jx) * p.*q;

nominal_p = Jpdot / Jx;
nominal_q = Jqdot / Jy;
nominal_r = Jrdot / Jz;


diff_p = observations(:,1) - nominal_p;
diff_q = observations(:,2) - nominal_q;
diff_r = observations(:,3) - nominal_r;

training_out_a = [diff_p diff_q diff_r];   % 1*3
