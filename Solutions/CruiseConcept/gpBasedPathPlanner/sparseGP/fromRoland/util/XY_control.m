function [ax_cmd, ay_cmd] = XY_control(X_REF, Y_REF, state, kp_xyz, kd_xyz)
    x_ref     = X_REF(1);
    vx_ref    = X_REF(2);
    ax_ref    = X_REF(3);        % feed forward
    x         = state(1);
    vx        = state(4);
    
    e_x       = x_ref - x;     % error
    e_vx      = vx_ref - vx;       % error derivative
    kp_x      = kp_xyz(1);           	% get kp and kd for x
    kd_x      = kd_xyz(1);
    ax_cmd    = kp_x*e_x + kd_x*e_vx + ax_ref;
    
    y_ref     = Y_REF(1);
    vy_ref    = Y_REF(2);
    ay_ref    = Y_REF(3);        % feed forward
    y         = state(2);
    vy        = state(5);
    
    e_y   = y_ref - y;     % error
    e_vy  = vy_ref - vy;   % error derivative
    kp_y  = kp_xyz(2);               % get kp and kd for y
    kd_y  = kd_xyz(2);
    ay_cmd = kp_y*e_y + kd_y*e_vy + ay_ref;
end