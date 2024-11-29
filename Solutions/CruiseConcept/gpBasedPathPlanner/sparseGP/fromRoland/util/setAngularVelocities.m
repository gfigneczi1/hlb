function [w1, w2, w3, w4] = setAngularVelocities(...
                                u1, pd_cmd, qd_cmd, rd_cmd, ...
                                state, I,para)
                            
    Ct = para.Ct;
    l = para.l;
    Cq = para.Cq;
    
    Ix = I(1);  Iy = I(2);  Iz = I(3);
    
    p = state(10);    q = state(11);    r = state(12);
    
    F_  = u1/ Ct;
    Mx_ = pd_cmd*Ix/(Ct*l) + (Iz-Iy)*q*r/(Ct*l);
    My_ = qd_cmd*Iy/(Ct*l) + (Ix-Iz)*p*r/(Ct*l);
    Mz_ = rd_cmd*Iz/Cq + (Iy-Ix)*p*q/Cq;
    
    w1 = (1/4)*(F_ + Mx_ + My_ + Mz_);
    w2 = (1/4)*(F_ - Mx_ + My_ - Mz_);
    w3 = (1/4)*(F_ - Mx_ - My_ + Mz_);
    w4 = (1/4)*(F_ + Mx_ - My_ - Mz_);
    
    w1 = -sqrt(w1); % since Tau1 (ccw) +ve
    w2 =  sqrt(w2);
    w3 = -sqrt(w3);
    w4 =  sqrt(w4);
end