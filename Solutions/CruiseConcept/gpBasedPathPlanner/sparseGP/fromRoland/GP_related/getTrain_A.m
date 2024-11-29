function training_in_a = getTrain_A(x,u,para)


% x,y,z, phi, theta, psi, vx, vy, vz, p,  q,  r
% 1 2 3   4     5     6   7   8   9   10  11  12

% state_hist-----  N*12
% F_hist --------  1*N


H = para.trans.H;


training_in_a = [x(10:12,1:end-H);u(4:6,:)]';  % 396*7