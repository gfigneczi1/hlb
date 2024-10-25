clear;
close all;

GENERATE_VIDEO = false;

TEST_CASES = [10 11 12 ...
    20 21 22 23 ...
    30 31 32 ...
    40 41 43 44 45 46 47];

%% Testcase definition
% 10 - constant low speed, no other objects
% 11 - constant mid speed, no other objects
% 12 - constant high speed, no other objects

% 20 - curvature based speed reduction, high speed in straights (normal deceleration), no other objects
% 21 - curvature based speed reduction, mid speed in straights (normal deceleration), no other objects
% 22 - curvature based speed reduction, mid speed in straights (intense deceleration), no other objects
% 23 - curvature based speed reduction, high speed in straights (intense deceleration), no other objects

% 30 - 10 with followed object truck - which decelerates in curve - obj:
% vxLow, eg: vxMid
% 31 - 11 with followed object small vehicle - which decelerates in curve -
% obj: vxMid, ego: vxHigh
% 32 - 11 with followed object small vehicle - which decelerates in curve -
% obj: vxHigh, ego: vxHigh

for t=1:length(TEST_CASES)
    TEST_CASE = TEST_CASES(t); 
    clear OncomingObjectType OncomingObjectPosition OncomingObjectDisplacement

% Test parameters definition
vxLow = 70/3.6; % m/s but given in km/h for more intuitive settings
vxMid = 90/3.6; % m/s but given in km/h for more intuitive settings
vxHigh = 110/3.6; % m/s but given in km/h for more intuitive settings

ayLow = 1; % m/s^2
ayMid = 2.5; % m/s^2 
ayHigh = 4; % m/s^2
axLow = 1; % in m/s^2 - acceleration
axMid = 2; % in m/s^2 - acceleration
axHigh = 3; % in m/s^2 - acceleration
daxLow = -1; % in m/s^2 - deceleration
daxMid = -2; % in m/s^2 - deceleration
daxHigh = -4; % in m/s^2 - deceleration

ObjectS0 = 20;

route = load("route.mat"); route = route.p;

%% static inputs: road curvature and curvature gradient
LaneOrientation = atan(movmean(diff(route(:,2))./diff(route(:,1)),300)); LaneOrientation = [LaneOrientation; LaneOrientation(end)]; % in radians
LaneCurvature = movmean(diff(tan(LaneOrientation))./diff(route(:,1)),200); LaneCurvature = [LaneCurvature; LaneCurvature(end)];
LaneCurvatureGradient = movmean(diff(tan(LaneCurvature))./diff(route(:,1)),200); LaneCurvatureGradient = [LaneCurvatureGradient; LaneCurvatureGradient(end)];

% calculate opposite lane
lw = 3.75;
routeOpp = [route(:,1)-sin(LaneOrientation)*lw route(:,2)+cos(LaneOrientation)*lw];
routeOpp = routeOpp(end:-1:1, :);
LaneOrientationOpp = atan(movmean(diff(routeOpp(:,2))./diff(routeOpp(:,1)),300)); LaneOrientationOpp = [LaneOrientationOpp; LaneOrientationOpp(end)]; % in radians
LaneCurvatureOpp = movmean(diff(tan(LaneOrientation))./diff(routeOpp(:,1)),200); LaneCurvatureOpp = [LaneCurvatureOpp; LaneCurvatureOpp(end)];
LaneCurvatureGradientOpp = movmean(diff(tan(LaneCurvatureOpp))./diff(routeOpp(:,1)),200); LaneCurvatureGradientOpp = [LaneCurvatureGradientOpp; LaneCurvatureGradientOpp(end)];
displacementOpp = cumtrapz(sqrt(diff(routeOpp(:,1)).^2+diff(routeOpp(:,2)).^2));
displacementOpp = [displacementOpp; displacementOpp(end)+displacementOpp(end)-displacementOpp(end-1)];

% refining the steps
%LaneOrientation = spline(route(:,1), LaneOrientation, 0:0.05:route(end,1))';
%LaneCurvature = spline(route(:,1), LaneCurvature, 0:0.05:route(end,1))';
route_y = spline(route(:,1), route(:,2), 0:0.05:route(end,1));
%route = [0:0.05:route(end,1); route_y]';
displacement = cumtrapz(sqrt(diff(route(:,1)).^2+diff(route(:,2)).^2));
displacement = [displacement; displacement(end)+displacement(end)-displacement(end-1)];

%% dynamic inputs - ego kinematic state: test case dependent
switch(TEST_CASE)
    case 10
        VelocityProfile = ones(size(route,1),1)*vxLow;
        AccelerationProfile = zeros(size(route,1),1);
        YawRateProfile = VelocityProfile.*LaneCurvature;
        FollowedObjectType = zeros(size(route,1),1);
        OncomingObjectType = zeros(size(route,1),1);
        EgoPosition = route;
        RelativeTime = displacement./VelocityProfile;
    case 11
        VelocityProfile = ones(size(route,1),1)*vxMid;
        AccelerationProfile = zeros(size(route,1),1);
        YawRateProfile = VelocityProfile.*LaneCurvature;
        FollowedObjectType = zeros(size(route,1),1);
        OncomingObjectType = zeros(size(route,1),1);
        EgoPosition = route;
        RelativeTime = displacement./VelocityProfile;
    case 12
        VelocityProfile = ones(size(route,1),1)*vxHigh;
        AccelerationProfile = zeros(size(route,1),1);
        YawRateProfile = VelocityProfile.*LaneCurvature;
        FollowedObjectType = zeros(size(route,1),1);
        OncomingObjectType = zeros(size(route,1),1);
        EgoPosition = route;
        RelativeTime = displacement./VelocityProfile;
    case 20
        % curvature based speed reduction, high speed in straights (normal deceleration), no other objects
        [RelativeTime, VelocityProfile, AccelerationProfile, ayProfile] = generateAccelerationLimitedMotion (vxHigh, axHigh, daxHigh, ayHigh, displacement, route, LaneCurvature);
        YawRateProfile = VelocityProfile.*LaneCurvature;
        FollowedObjectType = zeros(size(route,1),1);
        OncomingObjectType = zeros(size(route,1),1);
        EgoPosition = route;
    case 21
        % 21 - curvature based speed reduction, mid speed in straights (normal deceleration), no other objects
        [RelativeTime, VelocityProfile, AccelerationProfile, ayProfile] = generateAccelerationLimitedMotion (vxMid, axMid, daxMid, ayMid, displacement, route, LaneCurvature);
        YawRateProfile = VelocityProfile.*LaneCurvature;
        FollowedObjectType = zeros(size(route,1),1);
        OncomingObjectType = zeros(size(route,1),1);
        EgoPosition = route;
    case 22
        [RelativeTime, VelocityProfile, AccelerationProfile, ayProfile] = generateAccelerationLimitedMotion (vxLow, axLow, daxLow, ayLow, displacement, route, LaneCurvature);
        YawRateProfile = VelocityProfile.*LaneCurvature;
        FollowedObjectType = zeros(size(route,1),1);
        OncomingObjectType = zeros(size(route,1),1);
        EgoPosition = route;
    case 23
        [RelativeTime, VelocityProfile, AccelerationProfile, ayProfile] = generateAccelerationLimitedMotion (vxHigh, axLow, daxLow, ayLow, displacement, route, LaneCurvature);
        YawRateProfile = VelocityProfile.*LaneCurvature;
        FollowedObjectType = zeros(size(route,1),1);
        OncomingObjectType = zeros(size(route,1),1);
        EgoPosition = route;
    case 30
        % generate object profile
        objType = 2;
        idx0 = find(displacement >= ObjectS0,1);
        ObjectPositionX = route(idx0:end,1);
        ObjectPositionY = route(idx0:end,2);
        [RelativeTime, ObjectVelocityProfile] = generateAccelerationLimitedMotion (vxLow, axLow, daxLow, ayLow, displacement(idx0:end), [ObjectPositionX ObjectPositionY], LaneCurvature(idx0:end));
        % object quantities return reduced sample size as the object is
        % shifted a little
        LaneCurvature = LaneCurvature(1:end-idx0+1);
        LaneCurvatureGradient = LaneCurvatureGradient(1:end-idx0+1);
        LaneOrientation = LaneOrientation(1:end-idx0+1,1);
        
        % generate follow object profile
        [VelocityProfile, AccelerationProfile, FollowedObjectType, EgoPosition] = followObject(ObjectVelocityProfile, RelativeTime, axLow, daxLow, vxMid, 10, ObjectPositionX, ObjectPositionY, route, displacement, objType);
        OncomingObjectType = zeros(length(RelativeTime),1);
        YawRateProfile = VelocityProfile.*LaneCurvature;
        if (GENERATE_VIDEO)
            [frames] = createVideoFrames (EgoPosition, [ObjectPositionX ObjectPositionY], route, routeOpp, VelocityProfile, ObjectVelocityProfile, zeros(size(EgoPosition,1),2), OncomingObjectType, FollowedObjectType);
            myVideo = VideoWriter(fullfile(pwd, strcat("Scenario_", num2str(TEST_CASE), "_video.avi")));
            myVideo.FrameRate = 1/mean(diff(RelativeTime));
            open(myVideo);
            writeVideo(myVideo,frames);
            close(myVideo);
        end
    case 31
        % generate object profile
        objType = 1;
        idx0 = find(displacement >= ObjectS0,1);
        ObjectPositionX = route(idx0:end,1);
        ObjectPositionY = route(idx0:end,2);
        [RelativeTime, ObjectVelocityProfile] = generateAccelerationLimitedMotion (vxMid, axLow, daxLow, ayLow, displacement(idx0:end), [ObjectPositionX ObjectPositionY], LaneCurvature(idx0:end));
        % object quantities return reduced sample size as the object is
        % shifted a little
        LaneCurvature = LaneCurvature(1:end-idx0+1);
        LaneCurvatureGradient = LaneCurvatureGradient(1:end-idx0+1);
        LaneOrientation = LaneOrientation(1:end-idx0+1,1);
        
        % generate follow object profile
        [VelocityProfile, AccelerationProfile, FollowedObjectType, EgoPosition] = followObject(ObjectVelocityProfile, RelativeTime, axLow, daxLow, vxHigh, 10, ObjectPositionX, ObjectPositionY, route, displacement, objType);
        OncomingObjectType = zeros(length(RelativeTime),1);
        YawRateProfile = VelocityProfile.*LaneCurvature;
        if (GENERATE_VIDEO)
            [frames] = createVideoFrames (EgoPosition, [ObjectPositionX ObjectPositionY], route, routeOpp, VelocityProfile, ObjectVelocityProfile, zeros(size(EgoPosition,1),2), OncomingObjectType, FollowedObjectType);
            myVideo = VideoWriter(fullfile(pwd, strcat("Scenario_", num2str(TEST_CASE), "_video.avi")));
            myVideo.FrameRate = 1/mean(diff(RelativeTime));
            open(myVideo);
            writeVideo(myVideo,frames);
            close(myVideo);
        end
    case 32
        % generate object profile
        objType = 1;
        idx0 = find(displacement >= ObjectS0,1);
        ObjectPositionX = route(idx0:end,1);
        ObjectPositionY = route(idx0:end,2);
        [RelativeTime, ObjectVelocityProfile] = generateAccelerationLimitedMotion (vxHigh, axLow, daxLow, ayLow, displacement(idx0:end), [ObjectPositionX ObjectPositionY], LaneCurvature(idx0:end));
        % object quantities return reduced sample size as the object is
        % shifted a little
        LaneCurvature = LaneCurvature(1:end-idx0+1);
        LaneCurvatureGradient = LaneCurvatureGradient(1:end-idx0+1);
        LaneOrientation = LaneOrientation(1:end-idx0+1,1);
        
        % generate follow object profile
        [VelocityProfile, AccelerationProfile, FollowedObjectType, EgoPosition] = followObject(ObjectVelocityProfile, RelativeTime, axLow, daxLow, vxHigh, 10, ObjectPositionX, ObjectPositionY, route, displacement, objType);
        OncomingObjectType = zeros(length(RelativeTime),1);
        YawRateProfile = VelocityProfile.*LaneCurvature;
        if (GENERATE_VIDEO)
            [frames] = createVideoFrames (EgoPosition, [ObjectPositionX ObjectPositionY], route, routeOpp, VelocityProfile, ObjectVelocityProfile, zeros(size(EgoPosition,1),2), OncomingObjectType, FollowedObjectType);
            myVideo = VideoWriter(fullfile(pwd, strcat("Scenario_", num2str(TEST_CASE), "_video.avi")));
            myVideo.FrameRate = 1/mean(diff(RelativeTime));
            open(myVideo);
            writeVideo(myVideo,frames);
            close(myVideo);
        end
    case 40
        VelocityProfile = ones(size(route,1),1)*vxLow;
        AccelerationProfile = zeros(size(route,1),1);
        YawRateProfile = VelocityProfile.*LaneCurvature;
        FollowedObjectType = zeros(size(route,1),1);
        RelativeTime = displacement./VelocityProfile;
        EgoPosition = route;
        [OncomingObjectType, OncomingObjectPosition, OncomingObjectDisplacement] = generateOncomingObject(EgoPosition, RelativeTime, vxMid, VelocityProfile, 0, 0, routeOpp, displacementOpp, 1);
        if (GENERATE_VIDEO)
            [frames] = createVideoFrames (EgoPosition, zeros(size(EgoPosition,1),2), route, routeOpp, VelocityProfile, zeros(size(EgoPosition,1),1), OncomingObjectPosition, OncomingObjectType, FollowedObjectType);
            myVideo = VideoWriter(fullfile(pwd, strcat("Scenario_", num2str(TEST_CASE), "_video.avi")));
            myVideo.FrameRate = 1/mean(diff(RelativeTime));
            open(myVideo);
            writeVideo(myVideo,frames);
            close(myVideo);
        end
    case 41
        % generate object profile
        objType = 2; % truck
        idx0 = find(displacement >= ObjectS0,1);
        ObjectPositionX = route(idx0:end,1);
        ObjectPositionY = route(idx0:end,2);
        [RelativeTime, ObjectVelocityProfile] = generateAccelerationLimitedMotion (vxLow, axLow, daxLow, ayLow, displacement(idx0:end), [ObjectPositionX ObjectPositionY], LaneCurvature(idx0:end));
        % object quantities return reduced sample size as the object is
        % shifted a little
        LaneCurvature = LaneCurvature(1:end-idx0+1);
        LaneCurvatureGradient = LaneCurvatureGradient(1:end-idx0+1);
        LaneOrientation = LaneOrientation(1:end-idx0+1,1);
        
        % generate follow object profile
        [VelocityProfile, AccelerationProfile, FollowedObjectType, EgoPosition] = followObject(ObjectVelocityProfile, RelativeTime, axLow, daxLow, vxMid, 10, ObjectPositionX, ObjectPositionY, route, displacement, objType);
        [OncomingObjectType, OncomingObjectPosition, OncomingObjectDisplacement] = generateOncomingObject(EgoPosition, RelativeTime, vxMid, VelocityProfile, 0, 0, routeOpp, displacementOpp, 1);

        YawRateProfile = VelocityProfile.*LaneCurvature;
        if (GENERATE_VIDEO)
            [frames] = createVideoFrames (EgoPosition, [ObjectPositionX ObjectPositionY], route, routeOpp, VelocityProfile, ObjectVelocityProfile, OncomingObjectPosition, OncomingObjectType, FollowedObjectType);
            myVideo = VideoWriter(fullfile(pwd, strcat("Scenario_", num2str(TEST_CASE), "_video.avi")));
            myVideo.FrameRate = 1/mean(diff(RelativeTime));
            open(myVideo);
            writeVideo(myVideo,frames);
            close(myVideo);
        end
    case 43
        % generate object profile
        VelocityProfile = ones(size(route,1),1)*vxLow;
        AccelerationProfile = zeros(size(route,1),1);
        YawRateProfile = VelocityProfile.*LaneCurvature;
        FollowedObjectType = zeros(size(route,1),1);
        RelativeTime = displacement./VelocityProfile;
        EgoPosition = route;

        [OncomingObjectType{1}, OncomingObjectPosition{1}, OncomingObjectDisplacement{1}] = generateOncomingObject(EgoPosition, RelativeTime, vxMid, VelocityProfile, 0, 0, routeOpp, displacementOpp, 1);
        [OncomingObjectType{2}, OncomingObjectPosition{2}, OncomingObjectDisplacement{2}] = generateOncomingObject(EgoPosition, RelativeTime, vxMid, VelocityProfile, 0, 8, routeOpp, displacementOpp, 1);
        [OncomingObjectType{3}, OncomingObjectPosition{3}, OncomingObjectDisplacement{3}] = generateOncomingObject(EgoPosition, RelativeTime, vxMid, VelocityProfile, 0, 16, routeOpp, displacementOpp, 1);

        OncomingObjectType = mergeMultipleObjectTypes(OncomingObjectType, RelativeTime);
        
        if (GENERATE_VIDEO)
            [frames] = createVideoFrames (EgoPosition, zeros(size(EgoPosition,1),2), route, routeOpp, VelocityProfile, zeros(size(EgoPosition,1),1), OncomingObjectPosition, OncomingObjectType, FollowedObjectType);
            myVideo = VideoWriter(fullfile(pwd, strcat("Scenario_", num2str(TEST_CASE), "_video.avi")));
            myVideo.FrameRate = 1/mean(diff(RelativeTime));
            open(myVideo);
            writeVideo(myVideo,frames);
            close(myVideo);
        end
     case 44
        % generate object profile
        VelocityProfile = ones(size(route,1),1)*vxLow;
        AccelerationProfile = zeros(size(route,1),1);
        YawRateProfile = VelocityProfile.*LaneCurvature;
        FollowedObjectType = zeros(size(route,1),1);
        RelativeTime = displacement./VelocityProfile;
        EgoPosition = route;

        [OncomingObjectType{1}, OncomingObjectPosition{1}, OncomingObjectDisplacement{1}] = generateOncomingObject(EgoPosition, RelativeTime, vxLow, VelocityProfile, 0, 0, routeOpp, displacementOpp, 2);
        [OncomingObjectType{2}, OncomingObjectPosition{2}, OncomingObjectDisplacement{2}] = generateOncomingObject(EgoPosition, RelativeTime, vxLow, VelocityProfile, 0, 1, routeOpp, displacementOpp, 1);
        [OncomingObjectType{3}, OncomingObjectPosition{3}, OncomingObjectDisplacement{3}] = generateOncomingObject(EgoPosition, RelativeTime, vxLow, VelocityProfile, 0, 2, routeOpp, displacementOpp, 2);

        OncomingObjectType = mergeMultipleObjectTypes(OncomingObjectType, RelativeTime);
        
        if (GENERATE_VIDEO)
            [frames] = createVideoFrames (EgoPosition, zeros(size(EgoPosition,1),2), route, routeOpp, VelocityProfile, zeros(size(EgoPosition,1),1), OncomingObjectPosition, OncomingObjectType, FollowedObjectType);
            myVideo = VideoWriter(fullfile(pwd, strcat("Scenario_", num2str(TEST_CASE), "_video.avi")));
            myVideo.FrameRate = 1/mean(diff(RelativeTime));
            open(myVideo);
            writeVideo(myVideo,frames);
            close(myVideo);
        end
    case 45
        objType = 1; % small vehicle
        idx0 = find(displacement >= ObjectS0,1);
        ObjectPositionX = route(idx0:end,1);
        ObjectPositionY = route(idx0:end,2);
        [RelativeTime, ObjectVelocityProfile] = generateAccelerationLimitedMotion (vxMid, axLow, daxLow, ayLow, displacement(idx0:end), [ObjectPositionX ObjectPositionY], LaneCurvature(idx0:end));
        % object quantities return reduced sample size as the object is
        % shifted a little
        LaneCurvature = LaneCurvature(1:end-idx0+1);
        LaneCurvatureGradient = LaneCurvatureGradient(1:end-idx0+1);
        LaneOrientation = LaneOrientation(1:end-idx0+1,1);
        
        % generate follow object profile
        [VelocityProfile, AccelerationProfile, FollowedObjectType, EgoPosition] = followObject(ObjectVelocityProfile, RelativeTime, axLow, daxLow, vxHigh, 10, ObjectPositionX, ObjectPositionY, route, displacement, objType);
        % generate object profile
        YawRateProfile = VelocityProfile.*LaneCurvature;

        [OncomingObjectType{1}, OncomingObjectPosition{1}, OncomingObjectDisplacement{1}] = generateOncomingObject(EgoPosition, RelativeTime, vxMid, VelocityProfile, 0, 0, routeOpp, displacementOpp, 2);
        [OncomingObjectType{2}, OncomingObjectPosition{2}, OncomingObjectDisplacement{2}] = generateOncomingObject(EgoPosition, RelativeTime, vxMid, VelocityProfile, 0, 1.5, routeOpp, displacementOpp, 1);
        [OncomingObjectType{3}, OncomingObjectPosition{3}, OncomingObjectDisplacement{3}] = generateOncomingObject(EgoPosition, RelativeTime, vxMid, VelocityProfile, 0, 3, routeOpp, displacementOpp, 2);

        OncomingObjectType = mergeMultipleObjectTypes(OncomingObjectType, RelativeTime);
        
        if (GENERATE_VIDEO)
            [frames] = createVideoFrames (EgoPosition, zeros(size(EgoPosition,1),2), route, routeOpp, VelocityProfile, zeros(size(EgoPosition,1),1), OncomingObjectPosition, OncomingObjectType, FollowedObjectType);
            myVideo = VideoWriter(fullfile(pwd, strcat("Scenario_", num2str(TEST_CASE), "_video.avi")));
            myVideo.FrameRate = 1/mean(diff(RelativeTime));
            open(myVideo);
            writeVideo(myVideo,frames);
            close(myVideo);
        end
    case 46
        % generate object profile
        VelocityProfile = ones(size(route,1),1)*vxLow;
        AccelerationProfile = zeros(size(route,1),1);
        YawRateProfile = VelocityProfile.*LaneCurvature;
        FollowedObjectType = zeros(size(route,1),1);
        RelativeTime = displacement./VelocityProfile;
        EgoPosition = route;
        
        for i=1:8
            [OncomingObjectType{i}, OncomingObjectPosition{i}, OncomingObjectDisplacement{i}] = generateOncomingObject(EgoPosition, RelativeTime, vxMid, VelocityProfile, 0, 0+(i-1)*2, routeOpp, displacementOpp, 1);
        end

        OncomingObjectType = mergeMultipleObjectTypes(OncomingObjectType, RelativeTime);
        
        if (GENERATE_VIDEO)
            [frames] = createVideoFrames (EgoPosition, zeros(size(EgoPosition,1),2), route, routeOpp, VelocityProfile, zeros(size(EgoPosition,1),1), OncomingObjectPosition, OncomingObjectType, FollowedObjectType);
            myVideo = VideoWriter(fullfile(pwd, strcat("Scenario_", num2str(TEST_CASE), "_video.avi")));
            myVideo.FrameRate = 1/mean(diff(RelativeTime));
            open(myVideo);
            writeVideo(myVideo,frames);
            close(myVideo);
        end
    case 47
        % generate object profile
        VelocityProfile = ones(size(route,1),1)*vxLow;
        AccelerationProfile = zeros(size(route,1),1);
        YawRateProfile = VelocityProfile.*LaneCurvature;
        FollowedObjectType = zeros(size(route,1),1);
        RelativeTime = displacement./VelocityProfile;
        EgoPosition = route;
        
        for i=1:8
            [OncomingObjectType{i}, OncomingObjectPosition{i}, OncomingObjectDisplacement{i}] = generateOncomingObject(EgoPosition, RelativeTime, vxMid, VelocityProfile, 300, 0+(i-1)*2, routeOpp, displacementOpp, 2);
        end

        OncomingObjectType = mergeMultipleObjectTypes(OncomingObjectType, RelativeTime);
        
        if (GENERATE_VIDEO)
            [frames] = createVideoFrames (EgoPosition, zeros(size(EgoPosition,1),2), route, routeOpp, VelocityProfile, zeros(size(EgoPosition,1),1), OncomingObjectPosition, OncomingObjectType, FollowedObjectType);
            myVideo = VideoWriter(fullfile(pwd, strcat("Scenario_", num2str(TEST_CASE), "_video.avi")));
            myVideo.FrameRate = 1/mean(diff(RelativeTime));
            open(myVideo);
            writeVideo(myVideo,frames);
            close(myVideo);
        end
end

input = [OncomingObjectType ...
        FollowedObjectType ...
        VelocityProfile ...
        AccelerationProfile ...
        YawRateProfile ...
        LaneCurvature ...
        LaneCurvatureGradient];

data.input = input;
data.egoPosition = EgoPosition;
data.laneOrientation = LaneOrientation;
data.relativeTime = RelativeTime;

variables = ["$o_t$", "$f_{ot}$", "$v_x$", "$a_x$", "$\omega$", "$\kappa$", "$d\kappa$"];
f = figure();
set(f,'defaulttextInterpreter','latex') ;
set(f, 'defaultAxesTickLabelInterpreter','latex');  
set(f, 'defaultLegendInterpreter','latex');
for i=1:size(input,2)
    subplot(size(input,2),1,i);
    plot(RelativeTime, input(:,i));
    ylabel(variables(i));
    if (i==size(input,2))
        xlabel("Time(s)");
    end
    grid on;
    set(gca,'FontSize', 14);
end

close(f);

inputArray{t} = input;
save(strcat("synthetic_", num2str(TEST_CASE), ".mat"), "data");
fprintf(strcat("Test case", {' '}, num2str(TEST_CASE)," finished\n"));
end

function [frames] = createVideoFrames (EgoPosition, ObjectPosition, route, routeOpp, VelocityProfile, ObjectVelocityProfile, OncomingObjectPosition, OncomingObjectType, FollowedObjectType)
    set(0,'DefaultFigureVisible','off');
    f = figure();
    set(f,'defaulttextInterpreter','latex') ;
    set(f, 'defaultAxesTickLabelInterpreter','latex');  
    set(f, 'defaultLegendInterpreter','latex');
    if (~all(ObjectPosition==0))
        ds = sqrt(sum((EgoPosition'-ObjectPosition').^2));
    else
        ds = zeros(size(EgoPosition,1),1);
    end
    if (iscell(OncomingObjectPosition))
        dso = sqrt(sum((EgoPosition'-OncomingObjectPosition{1}').^2));
    else
        if (~all(OncomingObjectPosition==0))
            dso = sqrt(sum((EgoPosition'-OncomingObjectPosition').^2));
        else
            dso = zeros(size(EgoPosition,1),1);
        end
    end
    for i=1:size(EgoPosition,1)
        hold off;
        plot (EgoPosition(i,1), EgoPosition(i,2), 'ko', 'MarkerSize',8, 'DisplayName','EgoVehicle');
        hold on;
        grid on;
        plot (route(:,1), route(:,2), 'k--', 'DisplayName','route');
        plot (routeOpp(:,1), routeOpp(:,2), 'k--', 'DisplayName','routeOpp');
        if (iscell(OncomingObjectPosition))
            for j=1:length(OncomingObjectPosition)
                plot (OncomingObjectPosition{j}(i,1), OncomingObjectPosition{j}(i,2), 'ro', 'MarkerSize',8, 'DisplayName','OncomingVehicle');
            end
        else
            if (~all(OncomingObjectPosition(:,1)==0))
                plot (OncomingObjectPosition(i,1), OncomingObjectPosition(i,2), 'ro', 'MarkerSize',8, 'DisplayName','OncomingVehicle');
            end
        end

        if (~all(ObjectPosition(:,1)==0))
            plot (ObjectPosition(i,1), ObjectPosition(i,2), 'ro', 'MarkerSize',8, 'DisplayName','FollowedVehicle');
        end
        
        if (FollowedObjectType(i) == 0)
            strFollow = "none";
        elseif (FollowedObjectType(i) == 1)
            strFollow = "small vehicle";
        elseif (FollowedObjectType(i) == 2)
            strFollow = "truck";
        end
        if (iscell(OncomingObjectType))
           if (OncomingObjectType{1}(i) > 0)
               strOncoming = "convoy";
           else
               strOncoming = "none";
           end
        elseif (OncomingObjectType(i) == 0)
            strOncoming = "none";
        elseif (OncomingObjectType(i) == 1)
            strOncoming = "small vehicle";
        elseif (OncomingObjectType(i) == 2)
            strOncoming = "truck";
        elseif (OncomingObjectType == 3)
            strOncoming = "convoy";
        end

        str = strcat("Oncoming object is", {' '}, strOncoming);
        str2 = strcat("Follow object is", {' '}, strFollow);
        text(100,400,str);
        text(100,370,str2);

        legend;
        xlim([min(route(:,1)), max(route(:,1))]);
        ylim([min(route(:,2)), max(route(:,2))]);
        xlabel("$X_{UTM}(m)$"); ylabel("$Y_{UTM}(m)$");
        title(strcat("$v_{x,ego}=", num2str(VelocityProfile(i)), "m/s$, $v_{x,obj}=", num2str(ObjectVelocityProfile(i)), "m/s$, $d_{x,f}=", num2str(ds(i)), "m$, $d_{x,o}=", num2str(dso(i)), "$"));
        axis equal;
        frames(i) = getframe(gcf);
    end

    close(f);
    set(0,'DefaultFigureVisible','on');
end


function [RelativeTime, VelocityProfile, axProfile, ayProfile] = generateAccelerationLimitedMotion (vxNominal, axMax, daxMax, ayMax, displacement, route, LaneCurvature)
     % steady state (constant velocity)
    VelocityProfile = ones(size(route,1),1)*vxNominal;
    % lateral acceleration limits
    N = 2;
    firstCycle = true;
    n = 1;
    ayProfile = zeros(size(route,1),1); axProfile = zeros(size(route,1),1);
    while((any(abs(ayProfile) > ayMax*1.01) | any(axProfile > axMax*1.01) | any(axProfile<daxMax*1.01) | firstCycle) && n < 100)
        clear VelocityProfileFiltered;
        firstCycle = false;
        ayProfile = VelocityProfile.^2.*LaneCurvature; % vx^2*kappa
        ayProfile = min(max(-ayMax, ayProfile), ayMax);
        VelocityProfile = sqrt(ayProfile./LaneCurvature);
        RelativeTime = diff(displacement)./VelocityProfile(1:end-1);
        RelativeTime = cumtrapz(RelativeTime); 
        RelativeTime = [RelativeTime; RelativeTime(end)+diff(RelativeTime(end-1:end))];
        VelocityProfileFiltered(1) = VelocityProfile(1);
        for i=2:length(VelocityProfile)
            % filtering
            dvx = VelocityProfile(i) - VelocityProfileFiltered(i-1);
            dt = RelativeTime(i)-RelativeTime(i-1);
            dax = dvx/dt;
            VelocityProfileFiltered(i,1) = VelocityProfileFiltered(i-1,1)+min(max(daxMax, dax), axMax)*dt;
        end
        axProfile = diff(VelocityProfileFiltered)./diff(RelativeTime); axProfile = [axProfile; axProfile(end)];
        ayProfile = VelocityProfileFiltered.^2.*LaneCurvature;
        VelocityProfile(1:end-N+1) = VelocityProfileFiltered(N:end);
        n = n+1;
    end
    axProfile = movmean(axProfile,10);
    ayProfile = movmean(ayProfile,10);
    if (n>1000)
        fprintf("Maximum number of iterations reached.\n");
    end
end

function [VelocityProfile, AccelerationProfile, FollowedObjectType, EgoPosition] = followObject(ObjectVelocityProfile, RelativeTime, axMax, daxMax, vxNominal, ds0, ObjectPositionX, ObjectPositionY, route, displacement, objType)
    % this function calculates the velocity profile so that followed object
    % would be considered
    VelocityProfile(1,1) = vxNominal;
    AccelerationProfile(1,1) = 0;
    EgoPosition(1,1:2) = [route(1,1) route(1,2)];
    EgoDisplacement = 0;
    for i=2:length(RelativeTime)
        dt = RelativeTime(i)-RelativeTime(i-1);
        dx = EgoPosition(i-1,1)-ObjectPositionX(i);
        if (dx > 0)
            % ego vehicle is ahead of object, no need to control
            VelocityProfile(i,1) = vxNominal;
            ax = 0;
            FollowedObjectType(i-1,1) = 0;
            fprintf("Collision occured!\n");
        else
            % followed object is ahead of ego vehicle
            % calculate distance and relative speed to followed vehicle
            ds = sqrt(sum((EgoPosition(i-1,1:2)-[ObjectPositionX(i) ObjectPositionY(i)]).^2));
            dv = ObjectVelocityProfile(i-1)-VelocityProfile(i-1);
            dds = ds-ds0;
            if (dds < 0)
                % we are inside the safety range, deceleration with maximum
                % deceleration
                ax = daxMax;
            else
                % we are still far from safety range, calculate
                % acceleration
                t = dds/(dv/2+VelocityProfile(i-1));
                ax = min(max(dv/t, daxMax),axMax); % limit acceleration
            end
            % calculate next velocity point
            VelocityProfile(i,1) = VelocityProfile(i-1)+ax*dt;
            if (ds < 2*ds0)
                FollowedObjectType(i-1,1) = objType;
            else
                FollowedObjectType(i-1,1) = 0;
            end
        end
        EgoDisplacement(i) = EgoDisplacement(i-1)+VelocityProfile(i)*dt;
        EgoPosition(i,1:2) =  [spline(displacement,route(:,1),EgoDisplacement(i)) spline(displacement,route(:,2),EgoDisplacement(i))];
        AccelerationProfile(i,1) = ax;
    end
    FollowedObjectType(end+1,1) = FollowedObjectType(end,1);
end

function [OncomingObjectType, OncomingObjectPosition, OncomingObjectDisplacement] = generateOncomingObject(EgoPosition, RelativeTime, vObj, VelocityProfile, ObjS0, ObjT0, routeOpp, displacementOpp, objType)
% this function generates an oncoming object
    idx0 = find(displacementOpp >= ObjS0,1);
    idxStepIn = find(RelativeTime>=ObjT0,1);
    idxStepIn = max(idxStepIn,2);
    idx0 = max(idx0,2);
    OncomingObjectPosition(1:idxStepIn-1,1) = routeOpp(idx0,1);
    OncomingObjectPosition(1:idxStepIn-1,2) = routeOpp(idx0,2);
    OncomingObjectDisplacement(1:idxStepIn-1) = ObjS0;
    for i=idxStepIn:length(RelativeTime)
        % object steps in the scenario
        dt = RelativeTime(i)-RelativeTime(i-1);        
        OncomingObjectDisplacement(i) = OncomingObjectDisplacement(i-1)+dt*vObj;
        if (OncomingObjectDisplacement(i) > max(displacementOpp))
            % end of road reached, freeze calculation
            OncomingObjectPosition(i,1:2) = OncomingObjectPosition(i-1,1:2);
        else
            OncomingObjectPosition(i,1:2) = [spline(displacementOpp, routeOpp(:,1), OncomingObjectDisplacement(i)) ...
                spline(displacementOpp, routeOpp(:,2), OncomingObjectDisplacement(i))];
        end
    end
    % Now calculate oncoming object type
    OncomingObjectType = zeros(length(OncomingObjectDisplacement),1);
    distanceToObj = sqrt(sum((OncomingObjectPosition'-EgoPosition').^2))';
    distanceLimit = (VelocityProfile+vObj)*2.5; % 2 seconds
    OncomingObjectType((distanceToObj <= distanceLimit) & OncomingObjectPosition(:,1)>=EgoPosition(:,1)) = objType;
end

function OncomingObjectTypeMerged = mergeMultipleObjectTypes(OncomingObjectType, RelativeTime)
    if (~iscell(OncomingObjectType))
        OncomingObjectTypeMerged = OncomingObjectType;
    else
        dt = mean(diff(RelativeTime)); % s
        Tahead = 1; % s
        for i=1:length(OncomingObjectType)
            ObjectArray(i,:) = OncomingObjectType{i}';
            ObjectArrayPulledForward(i,:) = ObjectArray(i,:);
            RisingEdges = find(diff(ObjectArray(i,:))>0);
            RisingEdgesPulledForward = max(1,RisingEdges-floor(Tahead/dt));
            ObjectArrayPulledForward(i,RisingEdgesPulledForward:RisingEdges) = ObjectArrayPulledForward(i,RisingEdges+1);
        end
        ObjectArrayMerged = sum(ObjectArrayPulledForward);
        startTraffic = 0;
        for i=2:length(ObjectArrayMerged)
            if (ObjectArrayMerged(i)>0 && ObjectArrayMerged(i-1)==0)
                % traffic appeared
                startTraffic = i;
            elseif (ObjectArrayMerged(i)==0 && startTraffic > 0)
                if (any(any(ObjectArrayPulledForward(:,startTraffic:i-1)==2)) && ...
                        any(any(ObjectArrayMerged(startTraffic:i-1) > 1)))
                    % there are multiple objects and at least one is a
                    % truck, set the traffic to convoy
                    ObjectArrayMerged(startTraffic:i-1) = 3;% convoy mixed
                else
                    ObjectArrayMerged(startTraffic:i-1) = 1;% convoy small vehicle
                end
                startTraffic = 0;
            end
        end
        OncomingObjectTypeMerged = ObjectArrayMerged';
    end
end