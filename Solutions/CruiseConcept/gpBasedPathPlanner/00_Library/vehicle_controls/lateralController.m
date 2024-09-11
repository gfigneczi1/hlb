function [targetSteeringAngle, lap] = lateralController(path, vehicleState, parameters)
% pure-pursuit
p_lookAheadTime = parameters.controller.lat;

lad = p_lookAheadTime*vehicleState.vx;
% convert path to local frame
T = [cos(vehicleState.theta) sin(vehicleState.theta); ...
    -sin(vehicleState.theta) cos(vehicleState.theta)];
localPath = (path - [vehicleState.X vehicleState.Y])*T';
localPath(localPath(:,1)<0,2) = 0;
localPath(localPath(:,1)<0,1) = 0;
% find look ahead point
idx = find((localPath(:,1).^2+localPath(:,2).^2).^0.5 >= lad,1);
L = localPath(idx(1,1),1);
y = spline(localPath(localPath(:,1)>0,1), localPath(localPath(:,1)>0,2), L);

lap = [L, y];

% pure-pursuit relation
targetSteeringAngle = atan(parameters.vehicleParameters.wheelBase*(2*y)/L^2);
end

