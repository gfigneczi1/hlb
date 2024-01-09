function [local_coordinates] = GPS2ego(ego_pose,GPS_coordinates)
ego_pose_vector = ones(size(GPS_coordinates,1),2) * [ego_pose(1) 0; 0 ego_pose(2)];

T = [cos(ego_pose(3)) sin(ego_pose(3)); -sin(ego_pose(3)) cos(ego_pose(3))];
local_coordinates = (T*(GPS_coordinates-ego_pose_vector)')';
end

