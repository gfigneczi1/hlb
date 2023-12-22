function animateLDM(segment, GT_U, GT_Y, traj, cor, ref, validPoints, replan)

%% prepare data for animation

%road and car data
time_ = segment.q_T0-segment.q_T0(1);

[corridor, center] = createAbsRoadLane(segment, 'c01_left', 'c01_right');
T = [cos(segment.theta_calc(1)) -sin(segment.theta_calc(1)); sin(segment.theta_calc(1)) cos(segment.theta_calc(1))];
segment.theta_calc = movmean(segment.theta_calc, 100);

carpos=[segment.X_abs-segment.X_abs(1) segment.Y_abs-segment.Y_abs(1)]*T;

widths = abs(segment.c01_right - segment.c01_left);
roadwidth = mean(widths(~isnan(widths)));

traj_ = (ref-[segment.X_abs(1) segment.Y_abs(1)])*T; 
cor_ = (cor-[segment.X_abs(1) segment.Y_abs(1)])*T; 

% for i = 1:length(segment.X_abs)
%     x = (0:0.1:segment.VelocityX_ESP(i)*1.5)'; %(0:0.1:output(1).simout.s_end_traj.Data(i))';
%     theta = segment.theta_calc(i)-segment.theta_calc(1);
%     T = [cos(theta) -sin(theta); sin(theta) cos(theta)]; 
% end



scenario = drivingScenario;
roadCenters = corridor;
solidW = laneMarking('Solid','Width',0.3);
dashW = laneMarking('Dashed','Space',5);
lspec = lanespec([1 1], 'Width', roadwidth, 'Marking', [solidW dashW solidW]);
road(scenario,roadCenters,'Lanes',lspec);
car = vehicle(scenario, 'ClassID', 1, 'Position', carpos(1));
carGhost = vehicle(scenario, 'ClassID', 2, 'Position', center(1));

chasePlot(car);

hold on

dt = 0.05;

% make a figure and reset it.
% hFigure = figure(1);
% clf(hFigure,'reset');
% 
% x0=10;
% y0=10;
% width=1280;
% height=780;
% set(gcf,'position',[x0,y0,width,height])
% 
% add a scenario plot 
% hAxes = subplot(1,2,1);
% plot(segment.X_abs, segment.Y_abs);
% hold on;
% plot(scenario,'Parent',hAxes)
% title('Scenario plot');
% 
% add an egocentric plot for the vehicle
% hPanel = uipanel(hFigure,'Position',[.5 0 .5 1],'Units','Normal');
% chasePlot(car,'Parent',axes(hPanel))
% title('Chase plot');
% 
% Start the simulation loop
% for i = 1:length(segment.c01_left)
%     car.Position = [carpos(i,:) 0];
%     car.Yaw = (segment.theta_calc(i)-segment.theta_calc(1))*180/pi;
%     carGhost.Position = [center(i,:) 0];
%     carGhost.Yaw = (segment.theta_calc(i)-segment.theta_calc(1))*180/pi;
%     add a scenario plot 
%     updatePlots(scenario);
%     trajLine = plot(traj_(:,1),traj_(:,2), 'g', 'LineWidth', 3);
%     plot(hAxes, segment.X_abs, segment.Y_abs, segment.X_abs(i), segment.Y_abs(i), 'bo');
%     
%     Frames(i) = getframe(gcf);    
% end



for i = 11000:13000
    car.Position = [carpos(i,:) 0];
    car.Yaw = (segment.theta_calc(i)-segment.theta_calc(1))*180/pi;
    carGhost.Position = [center(i,:) 0];
    carGhost.Yaw = (segment.theta_calc(i)-segment.theta_calc(1))*180/pi;
    
%     if i > 25 && car.Position(2) < 4
%         car.Position = car.Position + [0 vy*dt 0];
%     end

    updatePlots(scenario)
    if (validPoints(i) == 1)
        trajLine = plot(traj_(:,1),traj_(:,2), 'g', 'LineWidth', 3);
        corLine = plot(cor_(:,1),cor_(:,2), 'b--', 'LineWidth', 2);
    end
%     if i <= length(pp_c00)
%         trajLine = plot(Trajectory{i}(:,1),Trajectory{i}(:,2), 'g');
%         corridorLine = plot(corridorTRP{i}(:,1), corridorTRP{i}(:,2), 'r');
%     end
    Frames(i-10999) = getframe(gcf);
%     if i <= length(pp_c00)
%         delete(trajLine);
%         delete(corridorLine);
%     end
end

myVideo = VideoWriter(fullfile(pwd, 'LDM_video.avi'));
myVideo.FrameRate = 25;
open(myVideo);
writeVideo(myVideo,Frames);
close(myVideo);

f = figure();
j = 1;
for i=11000:13000
    if (replan(i) > 0 && size(replan,1) <= i)
        subplot(3,1,1);
        p = plot(GT_U(1,1:j), GT_Y(2,1:j), 'o', GT_U(1,j), GT_Y(2,j), 'bo');
        set(p(1), 'color', [204 229 255]/255);
        xlim([-4, 4]); ylim([-1,1]);
        grid on;
        
        subplot(3,1,2);
        p = plot(GT_U(1,1:j), GT_Y(3,1:j), 'o', GT_U(1,j), GT_Y(3,j), 'bo');
        set(p(1), 'color', [204 229 255]/255);
        xlim([-4, 4]); ylim([-1,1]);
        grid on;
        
        subplot(3,1,3);
        p = plot(GT_U(1,1:j), GT_Y(4,1:j), 'o', GT_U(1,j), GT_Y(4,j), 'bo');
        set(p(1), 'color', [204 229 255]/255);
        xlim([-4, 4]); ylim([-1,1]);
        grid on;
        j = j +1;
    end
    Frames(i) = getframe(gcf);
end
close(f);

myVideo = VideoWriter(fullfile(pwd, 'LDM_correlation_video.avi'));
myVideo.FrameRate = 25;
open(myVideo);
writeVideo(myVideo,Frames);
close(myVideo);

