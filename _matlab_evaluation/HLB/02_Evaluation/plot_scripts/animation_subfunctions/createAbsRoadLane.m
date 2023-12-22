function [corridor, center] = createAbsRoadLane(subsegment, leftlane, rightlane)
%CREATEABSROADLANE Summary of this function goes here
%   Detailed explanation goes here

X_current = subsegment.X_abs(1);
Y_current = subsegment.Y_abs(1);

% preparation of corridor, orientation and curvature 
for j=1:size(subsegment.X_abs,1)
    corridor(j,1:2) = pos_tf2GPS(subsegment.X_abs(j),subsegment.Y_abs(j),subsegment.theta_calc(j),(subsegment.(leftlane)(j)));
    center(j,1:2) = pos_tf2GPS(subsegment.X_abs(j),subsegment.Y_abs(j),subsegment.theta_calc(j),(subsegment.(leftlane)(j)+subsegment.(rightlane)(j))/2);
end
T = [cos(subsegment.theta_calc(1)) -sin(subsegment.theta_calc(1)); sin(subsegment.theta_calc(1)) cos(subsegment.theta_calc(1))];
corridor = (corridor-[X_current Y_current])*T;
center = (center-[X_current Y_current])*T;


end

