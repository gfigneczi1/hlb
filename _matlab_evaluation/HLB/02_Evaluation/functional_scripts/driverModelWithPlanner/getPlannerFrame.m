function [plannerFrame] = getPlannerFrame(segmentsParam, previousPlannerFrame, trajectoryPreviousPlannerFrame, replanCycle, anchorPoints, anchorIndeces)
% function claculates the plannerFrame coordinates in GPS coordinate system
% outputs: 
%   plannerFrame - the (X, Y, theta) of selected planner frame in GPS
%   (global) coordinate system.
% inputs:
%   segmentsParam - c_0, c_1 parameters of the planned curves
%   previousPlannerFrame - the (X, Y, theta) of previous planner frame
%   trajectoryPreviousPlannerFrame - X,Y coordinates of the previously planned trajectory
%   replanCycle - number of cycle to do replan (select new planning point)
%   anchorPoints - border point information of planned trajectories
%   anchorIndeces - indeces of border points
%
% @copyright (C) 2022 Robert Bosch GmbH.
% The reproduction, distribution and utilization of this file as
% well as the communication of its contents to others without express
% authorization is prohibited. Offenders will be held liable for the
% payment of damages. All rights reserved in the event of the grant
% of a patent, utility model or design.
% @version 0.1
% @date 2022-04-20

nbAnchorPoints = floor(size(anchorPoints,1)/3);
plannerFrame(1) = trajectoryPreviousPlannerFrame(replanCycle, 1);
plannerFrame(2) = trajectoryPreviousPlannerFrame(replanCycle, 2);

% find the relevant segment
segment_id = find( anchorPoints > plannerFrame(1), 1 ) - 1;
if ((segment_id==0) ||(segment_id > nbAnchorPoints))
    segment_id = nbAnchorPoints-1;
end

% get the parameters of selected segments
c0_selected = segmentsParam(segment_id, 4);
c1_selected = segmentsParam(segment_id, 5);
theta_0_selected = anchorPoints(2*nbAnchorPoints + segment_id);

% calculate arc length for Theta value of planner point
s_planner = arc_length(trajectoryPreviousPlannerFrame, anchorIndeces(segment_id), replanCycle);

% Theta calcuation for planner point
plannerFrame(3) = theta_0_selected + (c0_selected * s_planner) + ((c1_selected / 2) * s_planner * s_planner);

T_inv =  [cos(previousPlannerFrame(3)) -sin(previousPlannerFrame(3)) 0; sin(previousPlannerFrame(3)) cos(previousPlannerFrame(3)) 0; 0 0 1];
 
plannerFrame = T_inv * plannerFrame';
plannerFrame = plannerFrame + previousPlannerFrame;

end

function [s] = arc_length(previousTrajectory, startIndex, replanCycle)

    dx = diff(previousTrajectory(startIndex:replanCycle,1));
    dy = diff(previousTrajectory(startIndex:replanCycle,1));
    ds = sqrt(dx.*dx + dy.*dy);
    s = sum(ds)*0.7109;

end
