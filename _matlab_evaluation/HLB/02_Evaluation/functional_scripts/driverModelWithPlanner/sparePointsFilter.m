function [valid_indeces] = sparePointsFilter(traj)
% this function searches for spare points (sporadic points) along the
% trajectory and returns the valid indeces. Parameter: maximum distances
% between points and minimum length of one standalone subsegment.
% This happens based on the distance between the neighbouring points

% @copyright (C) 2022 Robert Bosch GmbH.
% The reproduction, distribution and utilization of this file as
% well as the communication of its contents to others without express
% authorization is prohibited. Offenders will be held liable for the
% payment of damages. All rights reserved in the event of the grant
% of a patent, utility model or design.
% @version 1.0
% @date 2022-12-07
distances = (diff(traj(:,1)).^2+diff(traj(:,2)).^2).^0.5;
valid_indeces = find(distances < 3);

valid_length = 10; % unit: m
critical_points = find((diff(traj(valid_indeces,1)).^2+diff(traj(valid_indeces,2)).^2).^0.5 > valid_length);
filtered_points = ones(length(valid_indeces),1);
for i=2:length(critical_points)
    if (critical_points(i) - critical_points(i-1)) < 100
        filtered_points(critical_points(i-1):critical_points(i)) = 0;
    end
end
valid_indeces = valid_indeces(filtered_points == 1);

end

