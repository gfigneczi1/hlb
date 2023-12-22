function [LCL_segment, LCR_segment] = invalid_lanechange_filter(segments)
%INVALID_LANECHANGE_FILTER Summary of this function goes here
%   Detailed explanation goes here
LCL_segment = segments.LCL_segments;
LCR_segment = segments.LCR_segments;
right_cut = double.empty;
left_cut = double.empty;
% if the diff between neighbour GPS time data bigger then 50 ms -> cut
% if the diff between neighbour intersected lane info data bigger
% then 0.3 m -> cut
for i=1:length(segments.LCL_segments)
    if ~isempty(find(diff(segments.LCL_segments(i).GPS_time) > 50, 1))
        left_cut = [left_cut i];
    elseif segments.LCL_segments(i).position(1) + segments.LCL_segments(i).position(end) == 3 && ...
            (~isempty(find(diff(segments.LCL_segments(i).c01_left) > 0.3, 1)) ...
            || ~isempty(find(diff(segments.LCL_segments(i).c01_left) < -0.3, 1)))
        left_cut = [left_cut i];
    elseif segments.LCL_segments(i).position(1) + segments.LCL_segments(i).position(end) < 3 && ...
            (~isempty(find(diff(segments.LCL_segments(i).c01_right) > 0.3, 1)) ...
            || ~isempty(find(diff(segments.LCL_segments(i).c01_right) < -0.3, 1)))
        left_cut = [left_cut i];
    elseif abs(segments.LCL_segments(i).c01_left(end) - segments.LCL_segments(i).c01_left(1)) < 3.38
        left_cut = [left_cut i];
    end
end

for i=1:length(segments.LCR_segments)
    if ~isempty(find(diff(segments.LCR_segments(i).GPS_time) > 50, 1))
        right_cut = [right_cut i];
    elseif segments.LCR_segments(i).position(1) + segments.LCR_segments(i).position(end) == 3 && ...
            (~isempty(find(diff(segments.LCR_segments(i).c01_left) > 0.3, 1)) ...
            || ~isempty(find(diff(segments.LCR_segments(i).c01_left) < -0.3, 1)))
        right_cut = [right_cut i];
    elseif segments.LCR_segments(i).position(1) + segments.LCR_segments(i).position(end) < 3 && ...
            (~isempty(find(diff(segments.LCR_segments(i).c01_right) > 0.3, 1)) ...
            || ~isempty(find(diff(segments.LCR_segments(i).c01_right) < -0.3, 1)))
        right_cut = [right_cut i];
    elseif abs(segments.LCR_segments(i).c01_left(end) - segments.LCR_segments(i).c01_left(1)) < 3.38
        right_cut = [right_cut i];
    end
end

if ~isempty(left_cut)
    LCL_segment(left_cut) = [];
end
if ~isempty(right_cut)
    LCR_segment(right_cut) = [];
end
end

