%% Initialize the LTI dynamics for quadcopter

%% Translational model
% x: x,y,z,vx,vy,vz;
% u: phi,theta,Ftotal
% static working point: phi=0; theta=0; psi=0; Ftotal=-mg;

LTI.trans.A = [ 0, 0, 0, 1, 0, 0;
                0, 0, 0, 0, 1, 0;
                0, 0, 0, 0, 0, 1;
                0, 0, 0, 0, 0, 0;
                0, 0, 0, 0, 0, 0;
                0, 0, 0, 0, 0, 0;];
            
LTI.trans.B =   [ 0,     0,         0;
                  0,     0,         0;
                  0,     0,         0;
                  0,  -para.g,      0;
                  para.g,  0,       0;
                  0,       0,   -1/para.m;];
              

% LTI.trans.B =   [ 0,     0,         0;
%                   0,     0,         0;
%                   0,     0,         0;
%                   1,     0,      0;
%                   0,     1,       0;
%                   0,       0,   1;];
              


LTI.trans.C = eye(6);
LTI.trans.D = zeros(6,3);

LTI_trans = ss(LTI.trans.A,LTI.trans.B,LTI.trans.C,LTI.trans.D);

LTI_trans  = c2d(LTI_trans,para.dT);

para.trans.P = dare(LTI_trans.A,LTI_trans.B,para.trans.Q,para.trans.R);
para.trans.P = 100*eye(6);


%% Rotational model
% x: phi,theta,psi,dphi,dtheta,dpsi
% u: U2,U3,U4
LTI.rot.A = [ 0, 0, 0, 1, 0, 0;
              0, 0, 0, 0, 1, 0;
              0, 0, 0, 0, 0, 1;
              0, 0, 0, 0, 0, 0;
              0, 0, 0, 0, 0, 0;
              0, 0, 0, 0, 0, 0;];
          
LTI.rot.B =   [    0,       0,     0;
                   0,       0,     0;
                   0,       0,     0;
                1/para.Jx,  0,     0;
                   0,  1/para.Jy,  0;
                   0,       0,  1/para.Jz ;];    
               
LTI.rot.C = eye(6);
LTI.rot.D = zeros(6,3);

LTI_rot = ss(LTI.rot.A,LTI.rot.B,LTI.rot.C,LTI.rot.D);
LTI_rot  = c2d(LTI_rot,para.dT);
para.rot.P = dare(LTI_rot.A,LTI_rot.B,para.rot.Q,para.rot.R);

             
               
