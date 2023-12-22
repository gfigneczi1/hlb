function [undisturbed_LC, enforced_LC] = undisturbed_laneChange_filter(LC_segment)
%UNDISTURBED_LANECHANGE_FILTER Summary of this function goes here
%   Detailed explanation goes here
undisturbed_LC = LC_segment;
enforced_LC = struct;

LongAcc_thd = 0.8;      % [m/s^2]
LongAcc_thd_break = -1;     % [m/s^2]
Velo_diff_thd = 17;     % [km/h]
cut_ind = double.empty;

for i=1:length(LC_segment)
    LongAcc = lowpass(LC_segment(i).AccelerationX_ESP,10,100);
    dv = abs(LC_segment(i).VelocityX_ESP(1)*3.6 -  ...
        LC_segment(i).VelocityX_ESP(end)*3.6); % velocity diff in km/h
    
    if (max(LongAcc) > LongAcc_thd && dv > Velo_diff_thd) || ...
            min(LongAcc) < LongAcc_thd_break
        cut_ind = [cut_ind i];
    end
end

if ~isempty(cut_ind)
    undisturbed_LC(cut_ind) = [];
    enforced_LC = LC_segment(cut_ind);
end
end

