function [segment, output, input] = functional_laneChangeOpenLoopResimulation(segment, type)
% This function replays the segment (recorded) scenario with VMC corridor
% prep and trajectory planner to see what kind of trajectory was planned.
addpath('C:\tsdct\ETAS\ASCET\V6_1_4\target\psl\simulink_block');
% addpath('C:\work\git\simulink_block\');

segment.X_abs = segment.LongPos_abs * 40075000 .* cos(segment.LatPos_abs*pi()/180) / 360;
segment.Y_abs = segment.LatPos_abs * 111.32*1000;

[input, c0, q_T0, blink_idx] = prepareInput(segment, type);
input(2:end,90) = input(1:end-1,90);
input_ = zeros(size(input,1),92);
input_(:,2:end) = input;
input_(:,1) = input_(:,2);
input_(isnan(input_)) = 0;
input_(isinf(input_)) = 0;
N = size(input,1);
Tsim = (N-1)*0.02;
assignin('base','input_',input_);
assignin('base','Tsim',Tsim);
out = sim('alc.slx');
output.simout = out.simout;

output.d0 = out.simout.d0.Data;
output.c0 = c0;
output.q_T0 = q_T0;
output.blink_idx = blink_idx;
disp('sim finished')
end

function [input, c0, q_T0, blink_idx] = prepareInput(segment,type)
% for this given PSL the following inputs must be provided:
% x_points_Left: 20 points in local CS for left lane edge
% y_points_Left: 20 points in local CS for left lane edge
% x_points_Right: 20 points in local CS for right lane edge
% y_points_Right: 20 points in local CS for right lane edge
% number of valid points (constant 20)
% yawRate (rad/sec)
% vxEgo (m/s)
% roadWheelAg (rad)
% BhvrID (=3 for ALC)
% ActivationRequest (=1 constant all along)
% BhvrStgy (=1 for Comfort)

% The inputs must be resampled to 20ms sample time and put into an array
% extended with the lane chanel starting from zero seconds
% This yields an N x 88 (1 time channel, 87 data channel) input array.
% N must be calculated as the following: N = lengthTimeOfManeuver / 0.02s +
% 1 (+1: time starting from zero).
%% 1. Calculating metadata
dt = mean(diff(segment.q_T0));
N = round((segment.q_T0(end)-segment.q_T0(1))/0.02)+1;
input = zeros(N,90);
if (std(diff(segment.q_T0)) <= 0.0051 && N<=length(segment.q_T0))
    q_T0_norm = segment.q_T0-segment.q_T0(1);
    input(:,1) = (0:0.02:(N-1)*0.02)';
    input(:,83) = interp1(q_T0_norm, movmean(segment.yawRateGPS*pi()/180,20),input(:,1));
    input(:,84) = interp1(q_T0_norm, segment.VelocityX_ESP,input(:,1)); 
    input(:,85) = interp1(q_T0_norm, segment.SteeringAngle,input(:,1)); 
    input(:,82) = 20*ones(N,1);
    input(:,86) = 3 * ones(N,1);
    input(:,87) = ones(N,1);
    input(:,88) = ones(N,1);
    
    % calculate lane edge points (suppose 90m look ahead distance)
    x = linspace(0,90,20);
    
    if strcmp(type, 'left')
        if segment.position(1) == 1
            c0_Left_target = interp1(q_T0_norm, segment.c02_left,input(:,1)); 
            c0_Right_target = interp1(q_T0_norm, segment.c01_left,input(:,1)); 
            c0_Left_orig = interp1(q_T0_norm, segment.c01_left,input(:,1)); 
            c0_Right_orig = interp1(q_T0_norm, segment.c01_right,input(:,1)); 
        elseif segment.position(1) == 0
            c0_Left_target = interp1(q_T0_norm, segment.c01_left,input(:,1)); 
            c0_Right_target = interp1(q_T0_norm, segment.c01_right,input(:,1)); 
            c0_Left_orig = interp1(q_T0_norm, segment.c01_right,input(:,1)); 
            c0_Right_orig = interp1(q_T0_norm, segment.c02_right,input(:,1)); 
        else 
            disp('bade starting position for left lane change')
        end
        blinkTime = q_T0_norm(find(segment.LCL_blinker, 1));
    elseif strcmp(type, 'right')
        if segment.position(1) == 1
            c0_Left_target = interp1(q_T0_norm, segment.c01_right,input(:,1)); 
            c0_Right_target = interp1(q_T0_norm, segment.c02_right,input(:,1)); 
            c0_Left_orig = interp1(q_T0_norm, segment.c01_left,input(:,1)); 
            c0_Right_orig = interp1(q_T0_norm, segment.c01_right,input(:,1)); 
        elseif segment.position(1) == 2
            c0_Left_target = interp1(q_T0_norm, segment.c01_left,input(:,1)); 
            c0_Right_target = interp1(q_T0_norm, segment.c01_right,input(:,1)); 
            c0_Left_orig = interp1(q_T0_norm, segment.c02_left,input(:,1)); 
            c0_Right_orig = interp1(q_T0_norm, segment.c01_left,input(:,1));
        else
            disp('bade starting position for right lane change')
        end
        blinkTime = q_T0_norm(find(segment.LCR_blinker, 1));
    else
        disp('bad type in function call asd...')
    end
    
    c1 = interp1(q_T0_norm, segment.c1,input(:,1)); 
    c2 = interp1(q_T0_norm, segment.c2,input(:,1)); 
    yLeft_target = c0_Left_target + c1*x + c2/2*x.^2;
    yRight_target = c0_Right_target + c1*x + c2/2*x.^2;
    yLeft_orig = c0_Left_orig + c1*x + c2/2*x.^2;
    yRight_orig = c0_Right_orig + c1*x + c2/2*x.^2;
    for j=1:N
        input(j,2:21) = x;
        input(j,42:61) = x;
    end
    input(:,22:41) = yLeft_target;
    input(:,62:81) = yRight_target;
    
    [~, blink_idx] = min(abs(input(:,1)-blinkTime));
    
    input(1:blink_idx,22:41) = yLeft_orig(1:blink_idx,:);
    input(1:blink_idx,62:81) = yRight_orig(1:blink_idx,:);
    
    %create output variables
    q_T0 = input(:,1);
    c0 = (c0_Left_target + c0_Right_target) / 2;
    
    k = 1;
    input(1,89) = 1;
    input(1,90) = 0;
    timeSinceLastUpdate = 0;
    if strcmp(type, 'right')
        crossIndex = 22;
    else
        crossIndex = 62;
    end
    for i = 2:length(input(:,crossIndex))
        timeSinceLastUpdate = timeSinceLastUpdate + 0.02;
        if abs(input(i,crossIndex) - input(i-1,crossIndex)) < 0.00125
            k = k + 1;
            if k > 10
                k = 1;
            end
            timeSinceLastUpdate = 0;
        end
        
        input(i,90) = timeSinceLastUpdate;
        input(i,89) = k;
    end
    input(:,91) = interp1(q_T0_norm, segment.AccelerationX_ESP,input(:,1)); 
end
end

