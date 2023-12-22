function [LCL, LCR, position] = lanechange_cutter(rawData)
%LANECHANGE_CUTTER Summary of this function goes here
%   Detailed explanation goes here
LCL = rawData.LCL_blinker;
LCR = rawData.LCR_blinker;
LC_act_L = find(LCL);
LC_act_R = find(LCR);

% defining params for maneuver seg
% TODO: move them to general param_struct
seg_params.lanechange_cutter.YawR_trh_start_left = 0.009*25; % [rad/s]
seg_params.lanechange_cutter.YawR_trh_start_right = 0.009*31; % [rad/s]
seg_params.lanechange_cutter.YawR_trh_mid = -0.0065*19; % [rad/s]
seg_params.lanechange_cutter.YawR_trh_end = -0.0054*36; % [rad/s]
seg_params.lanechange_cutter.settling_time = 140;  % [s] value * sampling rate (10 ms currently)
seg_params.lanechange_cutter.data_size = size(LCL);

% extending LC start with 0.8 sec, end with 3.5 sec
if LC_act_L(1) <= 80
    LCL(1:LC_act_L(1)) = 1;
else
    LCL(LC_act_L(1)-80:LC_act_L(1)) = 1;
end

if LC_act_R(1) <= 80
    LCR(1:LC_act_R(1)) = 1;
else
    LCR(LC_act_R(1)-80:LC_act_R(1)) = 1;
end

for i=2:length(LC_act_L)
    if LC_act_L(i) - LC_act_L(i-1) > 1
        LCL(LC_act_L(i-1):LC_act_L(i-1)+350) = 1;
        LCL(LC_act_L(i)-80:LC_act_L(i)) = 1;
    end
end
if LC_act_L(end) ~= length(rawData.LCL_blinker)
    if length(rawData.LCL_blinker) - LC_act_L(end) > 350
        ext = 350;
    elseif length(rawData.LCL_blinker) - LC_act_L(end) < 350
        ext = length(rawData.LCL_blinker) - LC_act_L(end) - 10;
    end
    LCL(LC_act_L(end):LC_act_L(end)+ext) = 1;
end
for i=2:length(LC_act_R)
    if LC_act_R(i) - LC_act_R(i-1) > 1
        LCR(LC_act_R(i-1):LC_act_R(i-1)+350) = 1;
        LCR(LC_act_R(i)-80:LC_act_R(i)) = 1;
    end
end
if LC_act_R(end) ~= length(rawData.LCR_blinker)
    if length(rawData.LCR_blinker) - LC_act_R(end) > 350
        ext = 350;
    elseif length(rawData.LCR_blinker) - LC_act_R(end) < 350
        ext = length(rawData.LCR_blinker) - LC_act_R(end) - 10;
    end
    LCR(LC_act_R(end):LC_act_R(end)+ext) = 1;
end

LC_ext_act_L = find(LCL);
LC_ext_act_R = find(LCR);

% defining function (corrected YawRate)to detect maneuver segment of LCs
y = (rawData.yawRateESP-rawData.c2.*rawData.VelocityX_ESP).*rawData.VelocityX_ESP;
filtered_y = lowpass(y, 1, 100); % 1 Hz pb filtering data sampled with 100 Hz

% filtering maneuver
% direction left is 1, right is -1
LCL_new = maneuver_segmentor(LC_ext_act_L, filtered_y, seg_params.lanechange_cutter, 1);
LCR_new = maneuver_segmentor(LC_ext_act_R, filtered_y, seg_params.lanechange_cutter, -1);

LCL = LCL_new; LCR = LCR_new;

% intersection handling
intersection = find(LCR & LCL);
if ~isempty(intersection)
    intersection_start = 1;
    for i=2:length(intersection)
        if intersection(i) - intersection(i-1) > 1
            len_intersection = length(intersection(intersection_start:(i-1)));
            if mod(len_intersection, 2) == 1
                if LCL(intersection(intersection_start) - 1) == 1
                    LCR(intersection(intersection_start:...
                        (fix((i-intersection_start)/2)))) = 0;
                    LCL(intersection((fix((i-intersection_start)/2):...
                        i-1))) = 0;
                elseif LCR(intersection(intersection_start) - 1) == 1
                    LCL(intersection(intersection_start:...
                        (fix((i-intersection_start)/2)))) = 0;
                    LCR(intersection((fix((i-intersection_start)/2):...
                        i-1))) = 0;
                end
            elseif mod(len_intersection, 2) == 0
                if LCL(intersection(intersection_start) - 1) == 1
                    LCR(intersection(intersection_start:...
                        ((i-intersection_start)/2))) = 0;
                    LCL(intersection(((i-intersection_start)/2:...
                        i-1))) = 0;
                elseif LCR(intersection(intersection_start) - 1) == 1
                    LCL(intersection(intersection_start:...
                        ((i-intersection_start)/2))) = 0;
                    LCR(intersection(((i-intersection_start)/2:...
                        i-1))) = 0;
                end
            end
            intersection_start = i;
        end
    end
    len_intersection = length(intersection(intersection_start:end));
    if mod(len_intersection, 2) == 1
        if LCL(intersection(intersection_start) - 1) == 1
            LCR(intersection(intersection_start:...
                (fix((i+1-intersection_start)/2)))) = 0;
            if fix((i+1-intersection_start)/2) == 0
                inters_ind = 1;
            else
                inters_ind = fix((i+1-intersection_start)/2);
            end
            LCL(intersection((inters_ind:...
                i))) = 0;
        elseif LCR(intersection(intersection_start) - 1) == 1
            LCL(intersection(intersection_start:...
                (fix((i+1-intersection_start)/2)))) = 0;
            LCR(intersection((fix((i+1-intersection_start)/2):...
                i))) = 0;
        end
    elseif mod(len_intersection, 2) == 0
        if LCL(intersection(intersection_start) - 1) == 1
            LCR(intersection(intersection_start:...
                ((i+1-intersection_start)/2))) = 0;
            LCL(intersection(((i+1-intersection_start)/2:...
                i))) = 0;
        elseif LCR(intersection(intersection_start) - 1) == 1
            LCL(intersection(intersection_start:...
                ((i+1-intersection_start)/2))) = 0;
            LCR(intersection(((i+1-intersection_start)/2:...
                i))) = 0;
        end
    end
end

% position calculation
position = zeros(seg_params.lanechange_cutter.data_size);
LC_edges = edge(LCL)*1; LC_edges(find(edge(LCR))) = 2;  % 1 for, 2 for right
LC_index = find(LC_edges);
position(1:LC_index(1)+fix(((LC_index(2)-LC_index(1))/2))) = 1;
id = 1;
for i=1:2:length(LC_index) - 3
    if LC_edges(LC_index(i)) == 1 % LC direction left
        id = id + 1;
        position(LC_index(i)+fix(((LC_index(i+1)-LC_index(i))/2)):...
            LC_index(i+3)-fix((LC_index(i+3)-LC_index(i+2))/2)) = id;
    elseif LC_edges(LC_index(i)) == 2 % LC direction right
        id = id - 1;
        position(LC_index(i)+fix(((LC_index(i+1)-LC_index(i))/2)):...
            LC_index(i+3)-fix((LC_index(i+3)-LC_index(i+2))/2)) = id;
    end
end
if LC_edges(LC_index(end)) == 1
    id = id + 1;
    position(LC_index(end)-fix((LC_index(end)-LC_index(end-1))/2):end) = id;
elseif LC_edges(LC_index(end)) == 2
    id = id - 1;
    position(LC_index(end)-fix((LC_index(end)-LC_index(end-1))/2):end) = id;
end

% position correction in case of the measerument didn't started from
% position 1
% TODO: examine data by number of lanes in the current road 
% (e.g M6 is a 2x2 road, therefore the only valid positions are 1 and 2)
if ~isempty(find(position == -1, 1)) && isempty((find(position > 1, 1)))
    position = position + 1;
elseif ~isempty(find(position == 3, 1)) && isempty((find(position < 1, 1)))
    position = position - 1;
elseif isempty(find(position > 1, 1)) && isempty(find(position < 0, 1))
    position = position + 1;
end

end

function [LC] = maneuver_segmentor(LC_ongoing, YawRate, segment_params, direction)
%MANEUVER_SEGMENTOR Summary of this function goes here
%   Detailed explanation goes here
LC = zeros(segment_params.data_size);
k = 1;
man_start_point = 0;
man_end_point = -600;
count = 0;
LC_edges = find(edge(LC_ongoing));
LC_edges = [LC_edges LC_ongoing(end)];
for i=1:length(LC_ongoing)
    if isempty(k)
        k = j;
    end
    if count > 0 && i < LC_edges(count)
        continue
    end
    man_settled = false;
    if direction == 1
        if YawRate(LC_ongoing(i)) > segment_params.YawR_trh_start_left*direction && (LC_ongoing(i)-man_end_point > 600)
            man_start_point = LC_ongoing(i);
            for j=i:length(LC_ongoing)
                if YawRate(LC_ongoing(j)) < segment_params.YawR_trh_mid*direction
                    for k=j+100:length(LC_ongoing)
                        if YawRate(LC_ongoing(k)) > segment_params.YawR_trh_end*direction
                            man_settled = true;
                        end
                        if man_settled
                            man_end_point = LC_ongoing(k);
                            break
                        end
                    end
                end
                if man_settled
                    if man_end_point+segment_params.settling_time > length(LC)
                        ext = length(LC) - man_end_point - 10;
                        LC(man_start_point:man_end_point+ext) = 1;
                    else
                        LC(man_start_point:man_end_point+segment_params.settling_time) = 1;
                    end
                    count = count + 1;
                    break
                end
            end
        end
    elseif direction == -1
        if YawRate(LC_ongoing(i)) < segment_params.YawR_trh_start_right*direction && (LC_ongoing(i)-man_end_point > 600)
            man_start_point = LC_ongoing(i);
            for j=i:length(LC_ongoing)
                if YawRate(LC_ongoing(j)) > segment_params.YawR_trh_mid*direction
                    for k=j+100:length(LC_ongoing)
                        if YawRate(LC_ongoing(k)) < segment_params.YawR_trh_end*direction
                            man_settled = true;
                        end
                        if man_settled
                            man_end_point = LC_ongoing(k);
                            break
                        end
                    end
                end
                if man_settled
                    if man_end_point+segment_params.settling_time > length(LC)
                        ext = length(LC) - man_end_point - 10;
                        LC(man_start_point:man_end_point+ext) = 1;
                    else
                        LC(man_start_point:man_end_point+segment_params.settling_time) = 1;
                    end
                    break
                end
            end
        end  
    end
end
end

