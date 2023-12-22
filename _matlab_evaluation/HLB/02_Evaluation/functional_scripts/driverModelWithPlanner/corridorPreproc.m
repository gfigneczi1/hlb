function [corridorPlannerFrame] = corridorPreproc(subsegment, plannerFrame)
global indexes
% function claculates the corridor coordinates in the current planner coordinate system
% outputs:
%   corridorPlannerFrame
%       corridorPlannerFrameX
%       corridorPlannerFrameY
%       corridorPlannerFrameOrientation
%       corridorPlannerFrameCurvature
% inputs:
%   subsegment - subsegment data (500m segment from the current ego position)
%   deltaPos - X, Y, Theta difference between planner frame and nearest ego
%   point
%
% @copyright (C) 2022 Robert Bosch GmbH.
% The reproduction, distribution and utilization of this file as
% well as the communication of its contents to others without express
% authorization is prohibited. Offenders will be held liable for the
% payment of damages. All rights reserved in the event of the grant
% of a patent, utility model or design.
% @version 0.1
% @date 2022-04-20

corridor = [subsegment(:,indexes.corrX) subsegment(:,indexes.corrY)];

T = [cos(plannerFrame(3)) -sin(plannerFrame(3)); sin(plannerFrame(3)) cos(plannerFrame(3))];
corridorPlannerFrame = (corridor - [plannerFrame(1) plannerFrame(2)])*T;
corridorPlannerFrame(:,3) = subsegment(:,indexes.orientationVector) - plannerFrame(3);
corridorPlannerFrame(:,4) = movmean(subsegment(:,indexes.curvatureVector), 100);

end