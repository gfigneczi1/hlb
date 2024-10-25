function [transformed_position] = pos_tf2GPS(X, Y, theta, c01_left)
%POS_LEFT_GPS Summary of this function goes here
%   Detailed explanation goes here
X_GPS = -sin(theta)*c01_left + X;
Y_GPS = cos(theta)*c01_left + Y;
transformed_position = [X_GPS, Y_GPS];
end

