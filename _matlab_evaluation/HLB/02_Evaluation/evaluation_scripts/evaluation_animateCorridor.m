function evaluation_animateCorridor(output, input, input_raw, borders)
%animateCorridor animates the corridor provided to refLat. It also animates
%the pseudo ground-truth refline, vehicle path and the compensated
%reference line provided by the reflat.
%
% Syntax
% ======
%   animateCorridor(output, input)
%
% Description
% ===========
%   -
%
% Input Argument
% ===============
%   'input': provided by the trinity OL test file, this is based on the
%   dataSet provided by the simenvironment.
%   'output': provided by the trinity OL test file, this is generated when
%   executing the test, and signals match the output interface of the OMCL
%
% Exemplary Function Call
% =======================
%   animateCorridor(myOutput, myInput);
%
%--------------------------------------------------------------------------
% Bosch EYM2
%--------------------------------------------------------------------------

% History
% =======
%   Written by igg2bp,         17-Oct-2022
%
%   The script structure and some of the functions are based on:
%   hop2abt,        21-Jun-2019 09:57:27
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Options read
options.sampleTime = 0.1;
options.startOfVideo = 0;
options.maximumLengthOfVideo = 18;
options.frameRate = 6;

    %% Input mapping , change only when interface names change!
    % ActManLat
    inpt.X_abs = input.X_abs.Data;
    inpt.Y_abs = input.Y_abs.Data;
    inpt.orientation = input.theta_calc.Data;
    % Metadata
    inpt.vidName = "video";
    t = input.X_abs.Time - input.X_abs.Time(1);


    %% reformat reference line data from measurement
    % use delay compensated corridors. Yellow is typically not included in
    % measurement
    inpt.staticBlueLeft.x = input_raw.LeftPoints_x.Data;
    inpt.staticBlueLeft.y = input_raw.LeftPoints_y.Data;
    inpt.staticBlueRight.x = input_raw.RightPoints_x.Data;
    inpt.staticBlueRight.y = input_raw.RightPoints_y.Data;
    Tscorridor = mean(diff(input_raw.LeftPoints_x.Time));
    
    inpt.borderPointsRight = borders.borders.Data(:,41);
    inpt.BlueBorderRight.x = borders.borders.Data(:,1:2:39);
    inpt.BlueBorderRight.y = borders.borders.Data(:,2:2:40);
    inpt.borderPointsLeft = borders.borders.Data(:,82);
    inpt.BlueBorderLeft.x = borders.borders.Data(:,42:2:80);
    inpt.BlueBorderLeft.y = borders.borders.Data(:,43:2:81);
    
    % adding goalpoint
    inpt.lap.x(:,1) = borders.lap.Data(1,1,:);
    inpt.lap.y(:,1) = borders.lap.Data(1,2,:);
    inpt.lap.t(:,1) = borders.lap.Data(1,3,:) * pi()/180;
    
    %% calculate planned trajectory from LDM
    inpt.trajectory.x(:,:) = output.trajectoryEgoFrame.Data(:,1,:);
    inpt.trajectory.y(:,:) = output.trajectoryEgoFrame.Data(:,2,:);
    
    %% calculate planned trajectroy from LDM in plannerframe
    inpt.trajectorypf.x(:,:) = output.trajectoryPlannerFrame.Data(:,1,:);
    inpt.trajectorypf.y(:,:) = output.trajectoryPlannerFrame.Data(:,2,:);
    
    %% calculate input corridor to LDM
    inpt.corridor.x(:,:) = output.corridorEgoFrame.corridorEgo_X.Data(:,1,:);
    inpt.corridor.y(:,:) = output.corridorEgoFrame.corridorEgo_Y.Data(:,1,:);
    
    %% calculate input corridor to LDM in plannerframe
    inpt.corridorpf.x(:,:) = output.corridorPlannerFrame.corridor_X.Data(:,1,:);
    inpt.corridorpf.y(:,:) = output.corridorPlannerFrame.corridor_Y.Data(:,1,:);
    
    %% calculate video sample time
    % Measurement sample time:
    Ts = mean(diff(t));
    % Video sample time
    if (rem(options.sampleTime,Ts) <=1e-8 && options.sampleTime >= Ts)
        TsVideo = options.sampleTime;
    else
        if (options.sampleTime < Ts)
            TsVideo = Ts;
            disp('Warning: too low video sample time, setting measurement sample time instead');
        else
            TsVideo = round(options.sampleTime/Ts)*Ts;
            disp('Warning: option sample time is not multiplication of measurement sampling time. Rounding is applied.');
            disp(strcat('Used video sample time:',{' '},num2str(TsVideo),'secs'));
        end
    end
    % converting to step size as index:
    stepIdx = TsVideo / Ts;
    
    %% calculate video start and end point
    startIdx = options.startOfVideo/Ts+1;
    if (startIdx < 1)
        startIdx = 1;
        disp('Improper option set for video start time. Using 0 sec instead.');
    elseif(startIdx > length(t))
        disp('Improper option set for video start time. Using 0 sec instead.');
    end
    endIdx = options.maximumLengthOfVideo/Ts;
    if (endIdx < 1 || endIdx > length(t))
        endIdx = length(t);
        disp('Improper option set for video length. Using maximum time instead.');
    end
    
    startIdx = int32(startIdx);
    endIdx = int32(endIdx);
    stepIdx = int32(stepIdx);
    %% animate
    if ~isempty(inpt.vidName)
        v = VideoWriter(fullfile('../../_temp',inpt.vidName),'MPEG-4');
        v.FrameRate = options.frameRate;
        open(v);
    end

    f = figure();
    set(gcf,'Position',[1 41 1920 1080]);
    set(gcf,'Color', 'white');
    set(f,'Visible','off');
    hold all;
    disp('Video generation ongoing...');
    tstart = tic;
    for k = startIdx:stepIdx:endIdx
        
        % Calculating local odometry of the ego vehicle
        odometry = calculateVehicleOdometryLocal(inpt, k);
        
        xmax = (max(inpt.staticBlueLeft.x(k,:)) + max(inpt.staticBlueRight.x(k,:))) / 2;

        % Calculating local odometry of the ego vehicle which is within the
        % look ahead distance
        odometry = calculateVehicleOdometryLad (odometry, xmax, k);

        xlabel('x_{inertial} [m]');
        ylabel('y_{inertial} [m]');

        hold on;

        % local odometry
        %plotLocalOdometry(odometry, ...
                %{'b:', 'linewidth', 1, 'displayname', 'Local odometry'});

        % corridor and reference
        %plotCorridorInInertial(inpt, k*(Ts/Tscorridor), {'-or', 'displayname', 'blue corridor'});
        %plotMidLaneInInertial(inpt, k*(Ts/Tscorridor), {'-xk', 'displayname', 'reference'});
        
        % corridor input
        plotCorridorInput(inpt, k, {'k', 'displayname', 'corridor input'});
        
        % trajectory look ahead point
        %plotGoalPoint(inpt, k*(Ts/Tscorridor), {'ko', 'displayname', 'look-ahead point', 'LineWidth', 3});
        
        % trajectory
        plotTrajectory(inpt, k, {'b', 'displayname', 'planned trajectory'});

        legend show;
        text(5,6,sprintf('time = %1.3f',t(k)));
        grid on;
        xlim([-20 110]);
        axis equal;
        set(gca, 'Ylim', [-5,5]);
        
        % clean up
        if ~isempty(inpt.vidName)
            drawnow;
            writeVideo(v,getframe(gcf));
        end
        clf;
        tstop = toc(tstart);
        if (k==startIdx)
            tEstimated = tstop*(endIdx-startIdx)/stepIdx;
            disp(strcat('Estimated length of video generation:',{' '},num2str(tEstimated),{' '},'secs'));
        end
    end

    if ~isempty(inpt.vidName)
        close(v);
    end
    close(f);
end

%% sub functions
function odometry = calculateVehicleOdometryLocal(inpt, k)
    odometry.X_local = (inpt.X_abs-inpt.X_abs(k)).*cos(inpt.orientation(k))+ ...
        (inpt.Y_abs-inpt.Y_abs(k)).*sin(inpt.orientation(k));
    odometry.Y_local = -(inpt.X_abs-inpt.X_abs(k)).*sin(inpt.orientation(k))+ ...
        (inpt.Y_abs-inpt.Y_abs(k)).*cos(inpt.orientation(k));
    odometry.orientation_local = inpt.orientation - inpt.orientation(k);
end

function [odometry] = calculateVehicleOdometryLad (odometry, xmax, k)
    lookAheadPointIdx = find(odometry.X_local>=xmax,1);
    if (isempty(lookAheadPointIdx))
        % there are caes when the refline package is shorter than the
        % hostvehicle movement states
        odometry.X_local_lad = odometry.X_local(k:end-1);
        odometry.Y_local_lad = odometry.Y_local(k:end-1);
        odometry.orientation_local_lad = odometry.orientation_local(k:end-1);
    else
        odometry.X_local_lad = odometry.X_local(k:lookAheadPointIdx(1,1));
        odometry.Y_local_lad = odometry.Y_local(k:lookAheadPointIdx(1,1));
        odometry.orientation_local_lad = odometry.orientation_local(k:lookAheadPointIdx(1,1));
    end
end

function plotCorridorInInertial(inpt, k, plotOpts)
    plot(inpt.staticBlueLeft.x(k,:), inpt.staticBlueLeft.y(k,:), plotOpts{:});
    plot(inpt.staticBlueRight.x(k,:), inpt.staticBlueRight.y(k,:), plotOpts{:});
end

function plotGoalPoint(inpt, k, plotOpts)
    plot(inpt.lap.x(k,1), inpt.lap.y(k,1), plotOpts{:});
    arrow_x = cos(inpt.lap.t(k,1))*linspace(0,1,5);
    arrow_y = sin(inpt.lap.t(k,1))*linspace(0,1,5);
    plot(arrow_x+inpt.lap.x(k,1), arrow_y+inpt.lap.y(k,1), 'k','LineWidth',0.5, 'displayname','look ahead orientation');
end

function plotTransformedCorridor(inpt, k, plotOpts)
    n = inpt.borderPointsLeft(k);
    plot(inpt.BlueBorderLeft.x(k,1:n), inpt.BlueBorderLeft.y(k,1:n), plotOpts{:});
    n = inpt.borderPointsRight(k);
    plot(inpt.BlueBorderRight.x(k,1:n), inpt.BlueBorderRight.y(k,1:n), plotOpts{:});
end

function plotCorridorInput(inpt, k, plotOpts)
    n = inpt.corridor.x(end,k);
    plot(inpt.corridor.x(1:n,k), inpt.corridor.y(1:n,k), plotOpts{:});
end

function plotTrajectory(inpt, k, plotOpts)
    n = inpt.trajectory.x(end,k);
    plot(inpt.trajectory.x(1:n,k), inpt.trajectory.y(1:n,k), plotOpts{:});
end

function plotMidLaneInInertial(inpt, k, plotOpts)
    plot(0.5*(inpt.staticBlueRight.x(k,:)+inpt.staticBlueLeft.x(k,:)), ...
        0.5*(inpt.staticBlueRight.y(k,:)+ inpt.staticBlueLeft.y(k,:)), plotOpts{:});
end

function plotLocalOdometry(odometry, plotOpts)
    plot(odometry.X_local_lad, odometry.Y_local_lad, plotOpts{:});
end