function animateLC(LCL, output)
%ANIMATELC Summary of this function goes here
%   Detailed explanation goes here


%% prepare data for animation

%road and car data
time_ = LCL(1).q_T0-LCL(1).q_T0(1);
LCL(1).X_abs = interp1(time_, LCL(1).X_abs, output(1).simout.pp_c00.Time);
LCL(1).Y_abs = interp1(time_, LCL(1).Y_abs, output(1).simout.pp_c00.Time);
LCL(1).theta_calc = interp1(time_, movmean(LCL(1).theta_calc,25), output(1).simout.pp_c00.Time);
LCL(1).c02_left = interp1(time_, LCL(1).c02_left, output(1).simout.pp_c00.Time);
LCL(1).c01_left = interp1(time_, LCL(1).c01_left, output(1).simout.pp_c00.Time);
LCL(1).quality = ~isnan(LCL(1).X_abs);
% applying structure filter
LCL(1) = structure_filter(LCL(1),"quality",1);
corridor = createAbsRoadLane(LCL(1), 'c02_left', 'c01_left');
T = [cos(LCL(1).theta_calc(1)) -sin(LCL(1).theta_calc(1)); sin(LCL(1).theta_calc(1)) cos(LCL(1).theta_calc(1))];

carpos=[LCL(1).X_abs-LCL(1).X_abs(1) LCL(1).Y_abs-LCL(1).Y_abs(1)]*T;

widths = abs(LCL(1).c02_left - LCL(1).c01_left);
roadwidth = mean(widths(~isnan(widths)));

%planned trajectory
pp_c00 = output(1).simout.pp_c00.Data;
pp_c01 = output(1).simout.pp_c01.Data;
pp_c02 = output(1).simout.pp_c02.Data/2;
pp_c03 = output(1).simout.pp_c03.Data/6;
pp_c04 = output(1).simout.pp_c04.Data/24;
pp_c05 = output(1).simout.pp_c05.Data/120;

Trajectory = {};
corridorTRP = {};
for i = 1:length(LCL(1).X_abs)
    x = (0:0.1:LCL.VelocityX_ESP(i)*1.5)'; %(0:0.1:output(1).simout.s_end_traj.Data(i))';
    plannedTraj = pp_c00(i) + pp_c01(i)*x + pp_c02(i)*x.^2 + pp_c03(i)*x.^3 + pp_c04(i)*x.^4 + pp_c05(i)*x.^5;
    theta = LCL(1).theta_calc(i)-LCL(1).theta_calc(1);
    T = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    plannedTraj_rotated = [x plannedTraj] * T';
    plannedTraj_shifted = plannedTraj_rotated + carpos(i,:);
    Trajectory{i} = plannedTraj_shifted;
    
    corridorTRP_temp = output.simout.c0.Data(i) + output.simout.c1.Data(i)*x + output.simout.c2.Data(i)/2*x.^2 + output.simout.c3.Data(i)/6 * x.^3;
    corridorTRP_temp = [x corridorTRP_temp] * T';
    corridorTRP_temp = corridorTRP_temp + carpos(i,:);
    corridorTRP{i} = corridorTRP_temp;
    
end



scenario = drivingScenario;
roadCenters = corridor;
solidW = laneMarking('Solid','Width',0.3);
dashW = laneMarking('Dashed','Space',5);
lspec = lanespec(1, 'Width', roadwidth, 'Marking', [solidW solidW]);
road(scenario,roadCenters,'Lanes',lspec);
car = vehicle(scenario, 'ClassID', 1, 'Position', carpos(1));
plot(scenario)
chasePlot(car);
hold on

dt = 0.04;

for i = 1:length(LCL(1).c01_left)
    car.Position = [carpos(i,:) 0];
    car.Yaw = LCL(1).theta_calc(i)*180/pi;
    
%     if i > 25 && car.Position(2) < 4
%         car.Position = car.Position + [0 vy*dt 0];
%     end

    updatePlots(scenario)
    if i <= length(pp_c00)
        trajLine = plot(Trajectory{i}(:,1),Trajectory{i}(:,2), 'g');
        corridorLine = plot(corridorTRP{i}(:,1), corridorTRP{i}(:,2), 'r');
    end
    Frames(i) = getframe(gcf);
    if i <= length(pp_c00)
        delete(trajLine);
        delete(corridorLine);
    end
end

myVideo = VideoWriter(fullfile(pwd, 'LC_video.avi'));
myVideo.FrameRate = 25;
open(myVideo);
writeVideo(myVideo,Frames);
close(myVideo);



end

