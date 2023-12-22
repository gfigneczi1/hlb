function [segments] = laneChangePositionCorrection(segments, LC_direction)

%threshold describes the maximum value difference in meters between two points 
%next to each other
threshold = 2;
max_shift_val = 4;
lane_change_direction = 0;
reference_border = "c01_left";
border_fields = ["c01_right","c02_right","c02_left"];

if LC_direction == 'LCR'
    lane_change_direction = -1;
else
    lane_change_direction = 1;
end

for i=1:length(segments)
    %if the difference between two elements is higher than threshold, we have to
    %reinitialize position values
    position_change_idx = find(abs(diff(segments(i).(reference_border)))>threshold,1) + 1;
    segments(i).position(1:position_change_idx-1) = segments(i).position(1);
    segments(i).position(position_change_idx:end) = segments(i).position(1) + lane_change_direction;
    
    for j=1:length(border_fields)
        %check right direction
        shift_val = (find(abs(diff(segments(i).(border_fields(j))(position_change_idx:position_change_idx+max_shift_val)))>threshold,1))*-1;
        if isempty(shift_val)
            %check left direction
            shift_val = find(flip(abs(diff(segments(i).(border_fields(j))(position_change_idx-max_shift_val:position_change_idx))))>threshold,1)-1;
        end  
        if ~isempty(shift_val)
            segments(i) = shiftLaneInfo(segments(i), border_fields(j), shift_val);
        end
    end
    segments(i) = lane_arrangement(segments(i));
end



