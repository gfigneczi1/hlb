clear; clc; close all;

load("C:\git\hlb\Solutions\CruiseConcept\gpBasedPathPlanner\sparseGP\Dr008_Dr027_input_GP.mat");

for i=1:length(segments.segments)
    segment = segments.segments(i).segment;
    name = segments.segments(i).name;
    driverId = str2num(name(3:5));

    %% measurable quantity selection
    % straights
    straights = abs(segment.LaneCurvature)<=2.5e-4;
    curves = abs(segment.LaneCurvature)>2.5e-4;
    freeDriving = segment.FrontTrafficType==0;
    % straight line offset
    c0 = -0.5*(segment.LaneEdgePositionLeft+segment.LaneEdgePositionRight);
    c0Straight = c0(straights);
    KPI(driverId).meanc0Straights = mean(c0Straight);
    KPI(driverId).medianc0Straights = median(c0Straight);
    KPI(driverId).stdc0Straights = std(c0Straight);

    % straight line yawRate and acceleration (due to drift-and-compensate behavior)
    yawRateStraight = segment.YawRate(straights);
    ayStraight = segment.Acceleration_Y(straights);

    % curve yawRate and acceleration (due to drift-and-compensate behavior)
    yawRateCurve = segment.YawRate(curves);
    ayCurve = segment.Acceleration_Y(curves);

    KPI(driverId).meanYawRateStraights = mean(abs(yawRateStraight));
    KPI(driverId).medianYawRateStraights = median(abs(yawRateStraight));
    KPI(driverId).stdYawRateStraights = std(abs(yawRateStraight));

    KPI(driverId).meanYawRateCurves = mean(abs(yawRateCurve));
    KPI(driverId).medianYawRateCurves = median(abs(yawRateCurve));
    KPI(driverId).stdYawRateCurves = std(abs(yawRateCurve));

    KPI(driverId).meanAyStraights = mean(abs(ayStraight));
    KPI(driverId).medianAyStraights = median(abs(ayStraight));
    KPI(driverId).stdAyStraights = std(abs(ayStraight));

    KPI(driverId).meanAyCurves = mean(abs(ayCurve));
    KPI(driverId).medianAyCurves = median(abs(ayCurve));
    KPI(driverId).stdAyCurves = std(abs(ayCurve));

    % speed preferences: overall, straights and curves (when there is no
    % vehicle in front!)
    speedFreeDriving = segment.VelocityX(freeDriving);
    speedFreeDrivingCurves = segment.VelocityX(freeDriving & curves);
    speedFreeDrivingStraights = segment.VelocityX(freeDriving & straights);

    KPI(driverId).meanSpeedCurves = mean(speedFreeDrivingCurves);
    KPI(driverId).medianSpeedCurves = median(speedFreeDrivingCurves);
    KPI(driverId).stdSpeedCurves = std(speedFreeDrivingCurves);

    KPI(driverId).meanSpeedStraight = mean(speedFreeDrivingStraights);
    KPI(driverId).medianSpeedStraight = median(speedFreeDrivingStraights);
    KPI(driverId).stdSpeedStraight = std(speedFreeDrivingStraights);

    %% generating plots
    f = figure(1);
    f.Position = [100 100 1050 750];
    set(f,'defaulttextInterpreter','latex') ;
    set(f, 'defaultAxesTickLabelInterpreter','latex');  
    set(f, 'defaultLegendInterpreter','latex');

    subplot(2,3,1); % straight line offset
    histogram(c0Straight, 'Normalization','probability');
    xlabel("$c_0(m)$");
    ylabel("Probability");
    title(strcat("Straight line offset dr.", {' '}, num2str(driverId)));
    set(gca,'FontSize', 12);
    xlim([-0.5, 0.5]);

    subplot(2,3,2); % lateral acceleration
    histogram(ayCurve, 'Normalization','probability', "DisplayName","Curves");
    hold on;
    histogram(ayStraight, 'Normalization','probability', 'DisplayName', "Straights");
    xlabel("$a_y(m/s^2)$");
    ylabel("Probability");
    title(strcat("Lateral acceleration dr.", {' '}, num2str(driverId)));
    legend;
    set(gca,'FontSize', 12);
    xlim([-3,3]);

    subplot(2,3,4); % speed preferences
    histogram(speedFreeDrivingCurves, 'Normalization','probability', "DisplayName","Curves");
    hold on;
    histogram(speedFreeDrivingStraights, 'Normalization','probability', 'DisplayName', "Straights");
    xlabel("$v_x(m/s)$");
    ylabel("Probability");
    title(strcat("Speed preference dr.", {' '}, num2str(driverId)));
    legend;
    set(gca,'FontSize', 12);
    xlim([18,30]);

    subplot(2,3,5); % yawrate preferences
    histogram(yawRateCurve, 'Normalization','probability', "DisplayName","Curves");
    hold on;
    histogram(yawRateStraight, 'Normalization','probability', 'DisplayName', "Straights");
    xlabel("$\omega_z(rad/s)$");
    ylabel("Probability");
    title(strcat("Yaw rate preference dr.", {' '}, num2str(driverId)));
    legend;
    set(gca,'FontSize', 12);
    xlim([-0.1, 0.1]);

    subplot(1,3,3)
    set(gca,'XTick', [])
    set(gca,'YTick', [])
    set(gca,'XLim', [0,1])
    set(gca,'YLim', [0,1])

    str = sprintf(append("Offsets straight line:\n mean: ", num2str(KPI(driverId).meanc0Straights), " m", ...
        "\n median: ", num2str(KPI(driverId).medianc0Straights), " m", ...
        "\n deviation: ", num2str(KPI(driverId).stdc0Straights), " m", ...
        "\nLateral acceleration curves:\n mean: ", num2str(KPI(driverId).meanAyCurves), " $m/s^2$", ...
        "\n median: ", num2str(KPI(driverId).medianAyCurves), " $m/s^2$", ...
        "\n deviation: ", num2str(KPI(driverId).stdAyCurves), " $m/s^2$", ...
        "\nLateral acceleration straights:\n mean: ", num2str(KPI(driverId).meanAyStraights), " $m/s^2$", ...
        "\n median: ", num2str(KPI(driverId).medianAyStraights), " $m/s^2$", ...
        "\n deviation: ", num2str(KPI(driverId).stdAyCurves), " $m/s^2$", ...
        "\nYawRate straights:\n mean: ", num2str(KPI(driverId).meanYawRateStraights), " $rad/s$", ...
        "\n median: ", num2str(KPI(driverId).medianYawRateStraights), " $rad/s$", ...
        "\n deviation: ", num2str(KPI(driverId).stdYawRateStraights), " $rad/s$", ...
        "\nYawRate curves:\n mean: ", num2str(KPI(driverId).meanYawRateCurves), " $rad/s$", ...
        "\n median: ", num2str(KPI(driverId).medianYawRateCurves), " $rad/s$", ...
        "\n deviation: ", num2str(KPI(driverId).stdYawRateCurves), " $rad/s$", ...
        "\nSpeed preferences straights:\n mean: ", num2str(KPI(driverId).meanSpeedStraight), " $m/s$", ...
        "\n median: ", num2str(KPI(driverId).medianSpeedStraight), " $m/s$", ...
        "\n deviation: ", num2str(KPI(driverId).stdSpeedStraight), " $m/s$", ...
        "\nSpeed preferences curves:\n mean: ", num2str(KPI(driverId).meanSpeedCurves), " $m/s$", ...
        "\n median: ", num2str(KPI(driverId).medianSpeedCurves), " $m/s$", ...
        "\n deviation: ", num2str(KPI(driverId).stdSpeedCurves), " $m/s$"));
    text(0.05,0.5, str, "FontSize",10);

    close(f);

end