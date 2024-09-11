function vehicleState = vehicleModel(vehicleState, dT, mode, parameters)
% THIS FUNCTION includes simplified single-track models of vehicle defined
% in the parameters struct.
% - in1: vehicleState: a scalar struct with relevant vehicle kinematic and dynamic
% data; in the first cycle, initialize vehicleState with preferred values,
% then in the following cycles just simply return the calculated
% vehicleState
% - in2: dT - sample time of the discrete calculations
% - in3: mode - "kinematic", "dynamic" or "dynamicSimplified"
% - in4: parameters: scalar struct of vehicle parameters

if(mode == "kinematic")    
    % kinematic bicycle model

    % calculating yaw-rate based on steering angle
    vehicleState.yawRate = tan(vehicleState.steeringAngle)/parameters.wheelBase * vehicleState.vx;
    vehicleState.theta = vehicleState.theta + dT*vehicleState.yawRate;
    % displacement
    dX = vehicleState.vx*dT*cos(vehicleState.theta);
    dY = vehicleState.vx*dT*sin(vehicleState.theta);
    % new vehicle position
    vehicleState.X = vehicleState.X+dX;
    vehicleState.Y = vehicleState.Y+dY;
    vehicleState.vy = 0;
    vehicleState.ax = 0;
    vehicleState.ay = vehicleState.vx*vehicleState.yawRate;
    
elseif (mode == "dynamic")    
    % dynamic bicycle model
    % only for high speeds
    % inputs:
    % - vehicleState.steeringAngle
    % - vehicleState.M_av/M_ah - breaking forces
    % - vehicleState.M_bv/M_bh - breaking forces
    % parameters:
    % r: radius of wheels in [m]
    % c_alfav/h: lateral slip coefficient [N/rad]
    % c_sv/h: longitudinal slip coefficient [N/%]
    % c_w: wind coefficient
    % rho_air: air density
    % A: preface
    % J: rotational intertia of the vehicle
    % m: mass of the vehicle
    % lv/lh: COG position from front and rear axle
    % Jwheel: rotiational interatia of the wheels
%     r = 0.309725; c_alfaf = 76500; c_sf = 250; c_alfar = 76500; c_sr = 250;
%     m = 1519;
%     Jwheel = 250; J = 1818;
%     A = 1.5; c_w = 0; rho_air = 1; lf = 1; lr = 1.5; 
    
%    vehicleState.s_v = -(vehicleState.v_vx_v - r*vehicleState.drho_v)/(max(abs(vehicleState.v_vx_v), r*vehicleState.drho_v)); % longitudinal slip front, wheel coordinate frame
    vehicleState.s_f = -(vehicleState.v_fx_v - parameters.r*vehicleState.drho_f)/(parameters.r*vehicleState.drho_f); % longitudinal slip front, wheel coordinate frame

%    vehicleState.alfa_v = -vehicleState.v_vy_v /(abs(r*vehicleState.drho_v)); % lateral slip front, wheel coordinate frame
    vehicleState.alfa_f = -vehicleState.v_fy_v /(parameters.r*vehicleState.drho_f); % lateral slip front, wheel coordinate frame

%    vehicleState.s_h = -(vehicleState.v_hx_v - r*vehicleState.drho_h)/(max(abs(vehicleState.v_hx_v), r*vehicleState.drho_h)); % longitudinal slip front, wheel coordinate frame
    vehicleState.s_r= -(vehicleState.v_rx_v - parameters.r*vehicleState.drho_r)/(parameters.r*vehicleState.drho_r); % longitudinal slip front, wheel coordinate frame

%    vehicleState.alfa_h = -vehicleState.v_hy_v /(abs(r*vehicleState.drho_h)); % lateral slip front, wheel coordinate frame
    vehicleState.alfa_r = -vehicleState.v_ry_v /(parameters.r*vehicleState.drho_r); % lateral slip front, wheel coordinate frame

    vehicleState.F_fx_v = parameters.c_sf*vehicleState.s_f; % longitudinal tyre force front, wheel coordinate frame
    vehicleState.F_fy_v = parameters.c_alfaf*vehicleState.alfa_f; % lateral tyre force front, wheel coordinate frame
    vehicleState.F_fx = cos(vehicleState.steeringAngle)*vehicleState.F_fx_v - sin(vehicleState.steeringAngle)*vehicleState.F_fy_v; % longitudinal tyre force front, vehicle frame
    vehicleState.F_fy = sin(vehicleState.steeringAngle)*vehicleState.F_fx_v + cos(vehicleState.steeringAngle)*vehicleState.F_fy_v; % lateral tyre force front, vehicle frame
    vehicleState.F_rx = parameters.c_sr*vehicleState.s_r; % longitudinal tyre force, rear wheel, vehicle coordinate frame == wheel coordinate frame
    vehicleState.F_ry = parameters.c_alfar*vehicleState.alfa_r; % lateral tyre force, rear wheel, vehicle coordinate frame == wheel coordinate frame
    vehicleState.F_wx = 0.5*parameters.c_w*parameters.rho_air*parameters.A*vehicleState.v_fx^2;
    vehicleState.F_wy = 0.5*parameters.c_w*parameters.rho_air*parameters.A*vehicleState.v_fy^2;
    
    % COG quantities    
    vehicleState.a_x = 1/parameters.m*(vehicleState.F_fx + vehicleState.F_rx - vehicleState.F_wx);
    vehicleState.a_y = 1/parameters.m*(vehicleState.F_fy + vehicleState.F_ry - vehicleState.F_wy);
    
    vehicleState.v_x = vehicleState.v_x + vehicleState.a_x*dT;
    vehicleState.v_y = vehicleState.v_y + vehicleState.a_y*dT;
    
    vehicleState.eps_sigma = 1/parameters.J*(parameters.lf*vehicleState.F_fy - parameters.lr*vehicleState.F_ry);
    vehicleState.yawRate = vehicleState.yawRate + dT*vehicleState.eps_sigma;
    
    % Baselink quantities
    vehicleState.ax = vehicleState.a_x;
    vehicleState.vx = vehicleState.v_x;
    vehicleState.ay = vehicleState.yawRate*vehicleState.vx;
    vehicleState.vy = 0;
    
    % transforming to front and rear wheel
    vehicleState.v_fx = vehicleState.v_x;
    vehicleState.v_fy = vehicleState.v_y + vehicleState.yawRate*parameters.lf;
    vehicleState.v_rx = vehicleState.v_x;
    vehicleState.v_ry = vehicleState.v_y - vehicleState.yawRate*parameters.lr;
    
    vehicleState.v_fx_v = cos(vehicleState.steeringAngle)*vehicleState.v_fx + sin(vehicleState.steeringAngle)*vehicleState.v_fy;
    vehicleState.v_fy_v = -sin(vehicleState.steeringAngle)*vehicleState.v_fx + cos(vehicleState.steeringAngle)*vehicleState.v_fy;
    vehicleState.v_rx_v = vehicleState.v_rx;
    vehicleState.v_ry_v = vehicleState.v_ry;
    
    vehicleState.ddrho_f = 1/parameters.Jwheel * (vehicleState.M_af - vehicleState.M_bf)*sign(vehicleState.drho_f) - parameters.r*vehicleState.F_fx_v;
    vehicleState.ddrho_r = 1/parameters.Jwheel * (vehicleState.M_ar - vehicleState.M_br)*sign(vehicleState.drho_r) - parameters.r*vehicleState.F_rx;
    
    vehicleState.drho_f = vehicleState.drho_f + vehicleState.ddrho_f*dT;
    vehicleState.drho_r = vehicleState.drho_r + vehicleState.ddrho_r*dT;
    
    % absolute frame quantities
    vehicleState.theta = vehicleState.theta + dT*vehicleState.yawRate;
    % displacement
    dX = vehicleState.v_rx*dT*cos(vehicleState.theta)-sin(vehicleState.theta)*vehicleState.v_ry*dT;
    dY = vehicleState.v_rx*dT*sin(vehicleState.theta)+cos(vehicleState.theta)*vehicleState.v_ry*dT;
    % new vehicle position
    vehicleState.X = vehicleState.X+dX;
    vehicleState.Y = vehicleState.Y+dY;
    
    vehicleState.v = (vehicleState.v_x^2 + vehicleState.v_y^2)^0.5;  
    
elseif (mode == "dynamicSimplified")
%     r = 0.309725; c_alfaf = 76500; c_sf = 250; c_alfar = 76500; c_sr = 250;
%     m = 1519;
%     Jwheel = 250; J = 1818;
%     lf = 1; lr = 1.5;     
    
    df = vehicleState.steeringAngle;
    vehicleState.F_fx = cos(df)*(-parameters.c_sf*(cos(df)*vehicleState.v_x+sin(df)*vehicleState.v_y+sin(df)*vehicleState.yawRate*parameters.lf-parameters.r*vehicleState.drho_f)/(parameters.r*vehicleState.drho_f))+ ...
        sin(df)*parameters.c_alfaf*(-sin(df)*vehicleState.v_x+cos(df)*vehicleState.v_y+cos(df)*vehicleState.yawRate*parameters.lf)/(parameters.r*vehicleState.drho_f);
    vehicleState.F_fy = -sin(df)*parameters.c_sf*(cos(df)*vehicleState.v_x+sin(df)*vehicleState.v_y+sin(df)*vehicleState.yawRate*parameters.lf-parameters.r*vehicleState.drho_f)/(parameters.r*vehicleState.drho_f) - ...
        cos(df)*parameters.c_alfaf*(-sin(df)*vehicleState.v_x+cos(df)*vehicleState.v_y+cos(df)*vehicleState.yawRate*parameters.lf)/(parameters.r*vehicleState.drho_f);
    vehicleState.F_rx = parameters.c_sr*(-(vehicleState.v_x-parameters.r*vehicleState.drho_r)/(parameters.r*vehicleState.drho_f));
    vehicleState.F_ry = parameters.c_alfar*(-(vehicleState.v_y-vehicleState.yawRate*parameters.lr)/(parameters.r*vehicleState.drho_r));
    vehicleState.a_x = 1/parameters.m*(vehicleState.F_fx + vehicleState.F_rx);
    vehicleState.a_y = 1/parameters.m*(vehicleState.F_fy + vehicleState.F_ry);
    vehicleState.v_x = vehicleState.v_x + vehicleState.a_x*dT;
    vehicleState.v_y = vehicleState.v_y + vehicleState.a_y*dT;
    vehicleState.yawRate = vehicleState.yawRate + dT/parameters.J*parameters.lf*vehicleState.F_fy - dT*parameters.lr/parameters.J*vehicleState.F_ry;
    vehicleState.drho_f = vehicleState.drho_f+dT*parameters.r*parameters.c_sf*(cos(df)*vehicleState.v_x+sin(df)*vehicleState.v_y+sin(df)*vehicleState.yawRate*parameters.lf-parameters.r*vehicleState.drho_f)/(parameters.r*vehicleState.drho_f);
    vehicleState.drho_r = vehicleState.drho_r+parameters.r*dT*parameters.c_sr*(vehicleState.vx-parameters.r*vehicleState.drho_r)/(parameters.r*vehicleState.drho_r);
    
    vehicleState.theta = vehicleState.theta + vehicleState.yawRate*dT;
    % displacement
    v_hx = vehicleState.v_x;
    v_hy = vehicleState.v_y - parameters.lr*vehicleState.yawRate;
    dX = v_hx*dT*cos(vehicleState.theta)-sin(vehicleState.theta)*v_hy*dT;
    dY = v_hx*dT*sin(vehicleState.theta)+cos(vehicleState.theta)*v_hy*dT;
    % new vehicle position
    vehicleState.X = vehicleState.X+dX;
    vehicleState.Y = vehicleState.Y+dY;
    
end

vehicleState.time = vehicleState.time + dT;
end



