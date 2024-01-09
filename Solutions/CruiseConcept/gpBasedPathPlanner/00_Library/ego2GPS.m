function [GPS_pose] = ego2GPS(ego_pose,local_coordinates)

ego_pose_vector = ones(size(local_coordinates,1),2) * [ego_pose(1) 0; 0 ego_pose(2)];

T = [cos(ego_pose(3)) -sin(ego_pose(3)); sin(ego_pose(3)) cos(ego_pose(3))];
GPS_pose = (T*local_coordinates')' + ego_pose_vector;

end

