function [phi_dot, theta_dot, psi_dot] = getEulerDeriv(xrot)
    phi     = xrot(1);
    theta   = xrot(2);
    psi     = xrot(3);
    
    p = xrot(4);
    q = xrot(5);
    r = xrot(6);
    
    R = [1 sin(phi)*tan(theta) cos(phi)*tan(theta);
        0 cos(phi) -sin(phi);
        0 sin(phi)*sec(theta) cos(phi)*sec(theta)];
    
    euler_angle_dot= R*[p q r]';
    phi_dot = euler_angle_dot(1);
    theta_dot = euler_angle_dot(2);
    psi_dot = euler_angle_dot(3);
end