function state_dot = dynsys_Plin_Anonlin_wind(t,x,u,LTI,para,i)
%-- Now only position model is added with uncertainties--% 

global flag_wind  winds

% seperate the varables
    xtrans = x(1:6);
    xrot   = x(7:12);
    utrans = u(1:3);  % phi, theta, F
    urot   = u(4:6);
    
    phi = xrot(1); theta = xrot(2); psi = xrot(3);
    p   = xrot(4); q     = xrot(5); r   = xrot(6);
   
    
    % Compute orientation angle derivaties
    [phi_dot, theta_dot, psi_dot] = getEulerDeriv(xrot);
    
    % Compute body rate velocity derivatives
    [p_dot, q_dot, r_dot] = getBodyVelocityDot(urot, para, xrot);
    
    RotZYX = rotBodytoWorld(phi, theta, psi);
    
    dx_trans = LTI.trans.A*xtrans + LTI.trans.B*utrans;
    
     v       = x(4:6);
    
    dx_rot   = [phi_dot;
                theta_dot;
                psi_dot;
                p_dot;
                q_dot;
                r_dot];
    

    % Compute states' derivatives
    state_dot = [dx_trans;
                 dx_rot] ;% 12*1 column
             
    
     %% add wind disturbances to nominal model 
      % wind-related parameters
            a11 = 0.25*para.m;
%              a11 =0.3*para.m;
%             b12 = 0.6*para.m;
            b12 =0.6*para.m;
            c12 = para.e*a11;
            d11 = para.e^2*a11+para.e*2*sqrt(para.m*para.g)*0.1;
            d33 = a11*(para.l^2);

            zz = [0 0 1]';
            zzx = [0 -1  0;
                   1  0  0;
                   0  0  0;];
            proj = eye(3) - zz*zz'; % projection to XY plane

            A    = a11*proj;     
            B    = b12*zzx;
            C    = -c12*zzx;
            D    = d11*proj + d33*zz*zz';   % caution: here has included the root thrust; refer to 2017 IFAC paper by Kai.
     
            
      if flag_wind == 1
             
             wind = winds(i,:)';
             drag_force  = -(RotZYX*A*RotZYX')*(v-wind)  -  RotZYX*B*[p q r]';   % refer to the 2017 IFAC paper by Kai.

             drag_torque = -C*RotZYX'*(v-wind)-D*[p q r]';

%            state_dot = state_dot + [zeros(1,3),drag_force'/para.m, zeros(1,3),drag_torque(1)/para.Jx,drag_torque(2)/para.Jy,drag_torque(3)/para.Jz]';
             state_dot = state_dot + [zeros(1,3),drag_force'/para.m, zeros(1,6)]';
      end       
             
             
             
end