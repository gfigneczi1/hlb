function training_in_p = getTrain_P(x,u,para)


% x,y,z, phi, theta, psi, vx, vy, vz, p,  q,  r
% 1 2 3   4     5     6   7   8   9   10  11  12

% state_hist-----  N*12
% F_hist --------  1*N

% training_in_p = [state_hist(:,4:9) F_hist'];   % N*7 
H = para.trans.H;

% training_in_p = [res.x.rot(1:3,1:end-1) ; res.x.trans(4:6,1:end-1);res.u.trans(3,:)]';

training_in_p = [x(7:9,1:end-H);x(4:6,1:end-H);u(3,:)]';  % 396*7