function [p_dot, q_dot, r_dot] = getBodyVelocityDot(tau, para, xrot)
    tx = tau(1);    ty = tau(2);    tz = tau(3);
    Ix = para.Jx;      Iy = para.Jy;      Iz = para.Jz;
    
    p = xrot(4);
    q = xrot(5);
    r = xrot(6);
    
    p_ = tx - (Iz-Iy) * q * r;
    q_ = ty - (Ix-Iz) * p * r;
    r_ = tz - (Iy-Ix) * q * p;
    
    p_dot = p_/Ix;
    q_dot = q_/Iy;
    r_dot = r_/Iz;
end