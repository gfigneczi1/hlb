function state_dot = dynsys_Plin_Anonlin(t,x,u,LTI,para)
   
% seperate the varables
    xtrans = x(1:6);
    xrot   = x(7:12);
    utrans = u(1:3);  % phi, theta, F
    urot   = u(4:6);
    
   
    
    % Compute orientation angle derivaties
    [phi_dot, theta_dot, psi_dot] = getEulerDeriv(xrot);
    
    % Compute body rate velocity derivatives
    [p_dot, q_dot, r_dot] = getBodyVelocityDot(urot, para, xrot);
    
    
    dx_trans = LTI.trans.A*xtrans + LTI.trans.B*utrans;
    
    dx_rot   = [phi_dot;
                theta_dot;
                psi_dot;
                p_dot;
                q_dot;
                r_dot];
    
 
 
 

    % Compute states' derivatives
    state_dot = [dx_trans;
                 dx_rot] ;% 12*1 column
end