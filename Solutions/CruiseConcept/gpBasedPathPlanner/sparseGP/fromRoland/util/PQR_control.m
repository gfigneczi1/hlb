function [u2, u3, u4] = PQR_control(p_cmd, q_cmd, r_cmd, p, q, r, kp_pqr)
    ep = p_cmd - p;
    kp_p = kp_pqr(1);
    u2 = kp_p * ep;
    
    eq = q_cmd - q;
    kp_q = kp_pqr(2);
    u3 = kp_q * eq;
    
    er = r_cmd - r;
    kp_r = kp_pqr(3);
    u4 = kp_r * er;
end