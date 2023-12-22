function [ego_path, ego_orientation] = referenceCalculation(segment_X_abs, segment_Y_abs, segment_theta_calc, plannerFrame)
    T = [cos(plannerFrame(3)) -sin(plannerFrame(3)); sin(plannerFrame(3)) cos(plannerFrame(3))];
    ego_path = ([segment_X_abs segment_Y_abs] - [plannerFrame(1) plannerFrame(2)])*T;
    ego_orientation = segment_theta_calc - plannerFrame(3);
end

