function [lanevalidity, LCL, LCR] = lanevaliditycheck(segment)
%LANEVALIDITYCHECK Summary of this function goes here
%   Detailed explanation goes here
    LCL = zeros(length(segment.c01_left),1);
    LCR = zeros(length(segment.c01_left),1);
    k = 0;
    timeout = 1200;
    for i=1:length(segment.c01_left)
        interrupt_flag = false;
        lanechange_invalid = false;
        count = 0;
        if i > k
            if segment.Left_Index(i) == 1
                for j=i:length(segment.c01_left)
                    if count == timeout
                        LCL(i:i+count) = 0;
                        lanechange_invalid = true;
                        break
                    end
                    if (segment.LaneChange_Approved(j) == 1) && (interrupt_flag == false)...
                            && (segment.Left_Index(j) == 0)
                        for k=j:length(segment.c01_left) - 1
                            if segment.LaneChange_Approved(k+1) == 0
                                interrupt_flag = true;
                                break
                            end
                        end
                    elseif segment.Right_Index(j) == 1
                        lanechange_invalid = true;
                    end
                    if interrupt_flag == true
                        break
                    end
                    count = count + 1;
                end
                if lanechange_invalid == false && j ~= length(LCL)
                    LCL(i:j) = 1;
                elseif lanechange_invalid == true
                    LCL(i:j) = 0;
                end
            elseif segment.Right_Index(i) == 1
                for j=i:length(segment.c01_left)
                    if count == timeout
                        LCR(i:i+count) = 0;
                        lanechange_invalid = true;
                        break
                    end
                    if (segment.LaneChange_Approved(j) == 1) && (interrupt_flag == false)...
                            && (segment.Right_Index(j) == 0)
                        for k=j:length(segment.c01_left) - 1
                            if segment.LaneChange_Approved(k+1) == 0
                                interrupt_flag = true;
                                break
                            end
                        end
                    elseif segment.Left_Index(j) == 1
                      lanechange_invalid = true;
                    end
                    if interrupt_flag == true
                        break
                    end
                    count = count + 1;
                end
                if lanechange_invalid == false && j ~= length(LCR)
                    LCR(i:j) = 1;
                elseif lanechange_invalid == true
                    LCR(i:j) = 0;
                end
            end
        end
    end
    % Lane validity calculation
    dt = 0.05;
    debouncer = 0;
    laneinvalid = 1;
    lanevalidity = zeros(length(segment.c01_left),1);
    for i=2:length(segment.c01_left)-2
       if ((abs(segment.c01_left(i) - segment.c01_right(i)) < 1.875 ...
                && abs(segment.c01_left(i-1) - segment.c01_right(i-1)) < 1.875 ...
                && abs(segment.c01_left(i+1) - segment.c01_right(i+1)) < 1.875 ...
                && abs(segment.c01_left(i+2) - segment.c01_right(i+2)) < 1.875) ...
                || (segment.c01_left(i) == segment.c01_right(i) ...
                && segment.c01_left(i-1) == segment.c01_right(i-1)) ...
                || LCL(i) == 1 ...
                || LCR(i) == 1 ...
                || abs(segment.c01_left(i))<0.1 ...
                || abs(segment.c01_right(i))<0.1)
                % lane invalid
            laneinvalid = 1;
            debouncer = round(2/dt); %2 sec debounce timer
        else
            if (debouncer <=0)
                laneinvalid = 0;
            else
                debouncer = debouncer - 1;
            end
        end
        lanevalidity(i) = 1 - laneinvalid;
    end
    lanevalidity = lanevalidity';
    LCL = LCL'; LCR = LCR';
end

