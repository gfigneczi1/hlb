function [segment] = shiftLaneInfo(segment,dataField, numOfPositions)
%SHIFTLANEINFO Summary of this function goes here
%   Detailed explanation goes here

segment.(dataField) = circshift(segment.(dataField),numOfPositions);
if numOfPositions < 0
    segment.(dataField)(end+numOfPositions+1:end) = segment.(dataField)(end+numOfPositions);
elseif numOfPositions > 0
    segment.(dataField)(1:numOfPositions) = segment.(dataField)(numOfPositions+1);
end
      



