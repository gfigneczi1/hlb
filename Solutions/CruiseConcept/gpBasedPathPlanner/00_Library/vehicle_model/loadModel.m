function vehicleState = loadModel(vehicleState)
    b = 0; %0.0214;
    vehicleState.M_bf = vehicleState.v_fx * b;
    vehicleState.M_br = vehicleState.v_rx * b;
end

