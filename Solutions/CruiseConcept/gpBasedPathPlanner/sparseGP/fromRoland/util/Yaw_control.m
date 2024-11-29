function [r_cmd] = Yaw_control(psi_ref, psi, kp_roll_pitch_yaw)

    kp_yaw = kp_roll_pitch_yaw(3);
    e_psi = psi_ref - psi;
    r_cmd = kp_yaw * e_psi;
end