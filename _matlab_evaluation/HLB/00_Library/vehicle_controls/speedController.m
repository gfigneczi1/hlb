function vehicleState = speedController(vehicleState, targetSpeed)
    setSpeed = targetSpeed;
    P = 0.1;
    speedError = setSpeed - vehicleState.vx;
    vehicleState.M_av = 0; %speedError * P;
end

