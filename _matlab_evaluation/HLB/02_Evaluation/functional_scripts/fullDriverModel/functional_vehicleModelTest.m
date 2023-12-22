function [vehiclePath, vehicleStates] = functional_vehicleModelTest(path, targetSpeed)

modelMode = "dynamic";
parameters.vehicleParameters.wheelBase = 2.7;
parameters.vehicleParameters.r = 0.309725; 
parameters.vehicleParameters.c_alfaf = 8600; 
parameters.vehicleParameters.c_sf = 23000; 
parameters.vehicleParameters.c_alfar = 8600; 
parameters.vehicleParameters.c_sr = 23000;
parameters.vehicleParameters.m = 1519;
parameters.vehicleParameters.Jwheel = 250; 
parameters.vehicleParameters.J = 1818;
parameters.vehicleParameters.A = 1.5; 
parameters.vehicleParameters.c_w = 0; 
parameters.vehicleParameters.rho_air = 1; 
parameters.vehicleParameters.lf = 1; 
parameters.vehicleParameters.lr = 1.5; 

parameters.controller.lat = 4.5;

vehiclePath = [];
initPose = [0,-1,0];
vehicleState = initVehicleState(initPose, targetSpeed, modelMode);
dT = 0.002;
i = 1;

while(~scenarioFinishChecker(path, vehicleState, parameters))
    vehicleStates{i} = vehicleState;
    if (modelMode == "dynamic")
        vehicleState = loadModel(vehicleState);
        vehicleState = speedController(vehicleState, targetSpeed);
    end
    vehicleState.steeringAngle = lateralController(path, vehicleState, parameters);            
    
    % vehicle model
    vehicleState = vehicleModel(vehicleState, dT, modelMode, parameters.vehicleParameters);

    vehiclePath = [vehiclePath; [vehicleState.X vehicleState.Y]];
    i = i+1;    
end

end

function scenarioFinished = scenarioFinishChecker(path, vehicleState, parameters)
    p_lookAheadTime = parameters.controller.lat+0.1;

    lad = p_lookAheadTime*vehicleState.vx;
    % convert path to local frame
    T = [cos(vehicleState.theta) sin(vehicleState.theta); ...
    -sin(vehicleState.theta) cos(vehicleState.theta)];
    localPath = (path - [vehicleState.X vehicleState.Y])*T';
    localPath(localPath(:,1)<0,2) = 0;
    localPath(localPath(:,1)<0,1) = 0;
    % find look ahead point
    idx = find((localPath(:,1).^2+localPath(:,2).^2).^0.5 >= lad,1);
    if (isempty(idx))
        scenarioFinished = 1;
    else
        scenarioFinished = 0;
    end
end

function vehicleState = initVehicleState(initPose, initV0, mode)
if (mode=="kinematic")
    vehicleState.X = initPose(1);
    vehicleState.Y = initPose(2);
    vehicleState.theta = initPose(3);
    vehicleState.yawRate = 0;
    vehicleState.vx = initV0;
    vehicleState.vy = 0;
    vehicleState.ax = 0;
    vehicleState.ay = vehicleState.vx*vehicleState.yawRate;
    vehicleState.time = 0;    
    vehicleState.steeringAngle = 0;
elseif (mode =="dynamic")
    lf = 1; lr = 1.5; r = 0.309725;
    
    vehicleState.X = initPose(1);
    vehicleState.Y = initPose(2);
    vehicleState.theta = initPose(3);
    vehicleState.v_x = initV0;
    vehicleState.v_y = 0;
    vehicleState.time = 0;
    vehicleState.yawRate = 0;
    vehicleState.steeringAngle = 0;
    vehicleState.v_fx = vehicleState.v_x;
    vehicleState.v_fy = vehicleState.v_y + vehicleState.yawRate*lf;
    vehicleState.v_rx = vehicleState.v_x;
    vehicleState.v_ry = vehicleState.v_y - vehicleState.yawRate*lr;
    vehicleState.drho_f = vehicleState.v_fx/r;
    vehicleState.drho_r = vehicleState.v_rx/r;
    vehicleState.v_fx_v = cos(vehicleState.steeringAngle)*vehicleState.v_fx + sin(vehicleState.steeringAngle)*vehicleState.v_fy;
    vehicleState.v_fy_v = -sin(vehicleState.steeringAngle)*vehicleState.v_fx + cos(vehicleState.steeringAngle)*vehicleState.v_fy;
    vehicleState.v_rx_v = vehicleState.v_rx;
    vehicleState.v_ry_v = vehicleState.v_ry;
    
    vehicleState.vx = (vehicleState.v_x^2 + vehicleState.v_y^2)^0.5;
    vehicleState.ay = vehicleState.v_x*vehicleState.yawRate;
    
    vehicleState.M_af = 0; vehicleState.M_ar = 0;
    vehicleState.M_bf = 0; vehicleState.M_br = 0;
elseif (mode =="dynamicSimplified")
    lf = 1; lr = 1.5; r = 0.309725;
    
    vehicleState.X = initPose(1);
    vehicleState.Y = initPose(2);
    vehicleState.theta = initPose(3);
    vehicleState.v_x = initV0;
    vehicleState.v_y = 0;
    vehicleState.time = 0;
    vehicleState.yawRate = 0;
    vehicleState.steeringAngle = 0;
    vehicleState.v_fx = vehicleState.v_x;
    vehicleState.v_fy = vehicleState.v_y + vehicleState.yawRate*lf;
    vehicleState.v_rx = vehicleState.v_x;
    vehicleState.v_ry = vehicleState.v_y - vehicleState.yawRate*lr;
    vehicleState.drho_f = vehicleState.v_fx/r;
    vehicleState.drho_r = vehicleState.v_rx/r;
    vehicleState.v_fx_v = cos(vehicleState.steeringAngle)*vehicleState.v_fx + sin(vehicleState.steeringAngle)*vehicleState.v_fy;
    vehicleState.v_fy_v = -sin(vehicleState.steeringAngle)*vehicleState.v_fx + cos(vehicleState.steeringAngle)*vehicleState.v_fy;
    vehicleState.v_rx_v = vehicleState.v_rx;
    vehicleState.v_ry_v = vehicleState.v_ry;
    
    vehicleState.vx = (vehicleState.v_x^2 + vehicleState.v_y^2)^0.5;
    vehicleState.ay = vehicleState.v_x*vehicleState.yawRate;
    
    vehicleState.M_af = 0; vehicleState.M_ar = 0;
    vehicleState.M_bf = 0; vehicleState.M_br = 0;
end
end

