function [indeces, X] = nodePointModel(corridor, P)
% P is a parameter of 3, which gives the distances from between the points.
% Parameters are here limited to have a minimum distance from the end of
% the corridor and from each other. Therefore the first point cannnot be
% closer to the corridor end than the 2 x minimum distance, second point as
% the 1 x minimum distance, and the last point of course cannot be greater
% than corridor length. Once the parameter is calculated, the node points
% are discretized to the corridor steps.

% @copyright (C) 2022 Robert Bosch GmbH.
% The reproduction, distribution and utilization of this file as
% well as the communication of its contents to others without express
% authorization is prohibited. Offenders will be held liable for the
% payment of damages. All rights reserved in the event of the grant
% of a patent, utility model or design.
% @version 1.0
% @date 2022-12-07
N = length(P); % number of node points

lengths = (corridor(:,1).^2+corridor(:,2).^2).^0.5;

n = size(corridor,1);
minStepDist = 10; % in meters
indeces = zeros(N+1,1);
X = zeros(N+1,1);

X(1) = corridor(1,1);
indeces(1) = 1;
L(1) = 0;

for i=1:N
    % iterate through P distances
    % first element of L is always zero in the planner frame
    L(i+1) = min(lengths(end)-(N-i+1)*minStepDist, max(L(i)+minStepDist, L(i)+P(i)*250));
%     L(2) = min(lengths(end)-3*minStepDist,max(minStepDist,P(1)*250));
%     L(3) = min(lengths(end)-2*minStepDist,max(L(2)+minStepDist,L(2)+P(2)*250));
%     L(4) = min(lengths(end)-minStepDist,max(L(3)+minStepDist,L(3)+P(3)*250));
end

for i=1:N
    indeces(i+1) = find(lengths>=L(i+1),1);
end

% indeces(2) = find(lengths>=L(2),1);
% indeces(3) = find(lengths>=L(3),1);
% indeces(4) = min(n,find(lengths>=L(4),1));
X(2:end) = corridor(indeces(2:end),1);

end

