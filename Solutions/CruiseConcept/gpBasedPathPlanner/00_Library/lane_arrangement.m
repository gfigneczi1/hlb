function [rawData] = lane_arrangement(rawData)
%LANE_ARRANGEMENT Summary of this function goes here
%   Detailed explanation goes here
for i=1:length(rawData.c01_left)
    if (rawData.position(i) == 1)
        arrangedLanes = [ rawData.c02_left(i); rawData.c01_left(i); rawData.c01_right(i); rawData.c02_right(i)];
    elseif (rawData.position(i) == 0)
        % inner left, middle left, middle right, outer right
        arrangedLanes = [ NaN; rawData.c02_left(i); rawData.c01_left(i); rawData.c01_right(i)];
    elseif (rawData.position(i) == 2)
        arrangedLanes = [ rawData.c01_left(i); rawData.c01_right(i); rawData.c02_right(i); NaN];
    else
        arrangedLanes = [NaN; NaN; NaN; NaN];
    end
    rawData.c02_left(i) = arrangedLanes(1);
    rawData.c01_left(i) = arrangedLanes(2);
    rawData.c01_right(i) = arrangedLanes(3);
    rawData.c02_right(i) = arrangedLanes(4);
end
end

