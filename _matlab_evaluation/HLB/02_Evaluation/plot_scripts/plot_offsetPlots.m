function f = plot_offsetPlots(segment_m,indexes,offsets, driverID)

    relSegment = [0, 1200];
    dx = mean(diff(segment_m(:,indexes.X_abs)));
    X_abs = segment_m(:,indexes.X_abs);

    f = figure();
    f.Position = [10 10 1000 1500];
    set(f,'defaulttextInterpreter','latex') ;
    set(f, 'defaultAxesTickLabelInterpreter','latex');  
    set(f, 'defaultLegendInterpreter','latex');

    subplot(3,1,1);
    straightSections = abs(segment_m(:,indexes.LaneCurvature))<2.5e-4;
    straightBoundingBoxes = zeros(numel(find(diff(straightSections) < 0))+numel(find(diff(straightSections) > 0))+1,1);
    if (straightSections(1))
        % start with straight
        straightBoundingBoxes(2:2:end) = find(diff(straightSections) < 0);
        straightBoundingBoxes(1:2:end) = [1; find(diff(straightSections) > 0)];
    elseif (straightSections(end))
        % end with straight
        straightBoundingBoxes(2:2:end) = [find(diff(straightSections) < 0); size(segment_m,1)];
        straightBoundingBoxes(1:2:end) = [find(diff(straightSections) > 0)];
    else
        straightBoundingBoxes(1:2:end-1) = find(diff(straightSections) > 0);
        straightBoundingBoxes(2:2:end-1) = find(diff(straightSections) < 0);
    end

    grid on; hold on;
    for j=1:length(straightBoundingBoxes)/2
        if (j==1)
            fill([X_abs(straightBoundingBoxes(2*j-1)); X_abs(straightBoundingBoxes(2*j)); X_abs(straightBoundingBoxes(2*j)); X_abs(straightBoundingBoxes(2*j-1))], ...
                [-1;-1;1;1], 'g', 'FaceAlpha', 0.2, 'DisplayName', 'Straight sections', 'LineStyle','none');
            hold on;
        else
            fill([X_abs(straightBoundingBoxes(2*j-1)); X_abs(straightBoundingBoxes(2*j)); X_abs(straightBoundingBoxes(2*j)); X_abs(straightBoundingBoxes(2*j-1))], ...
                [-1;-1;1;1], 'g', 'FaceAlpha', 0.2, 'LineStyle','none', "HandleVisibility","off");
        end
    end
    for i=1:length(offsets)
        plot(offsets{i}.X, offsets{i}.offset, offsets{i}.marker,  'LineWidth', 1.5, 'DisplayName', offsets{i}.name);
    end
    xlabel('X-UTM(m)'); ylabel("$\delta(m)$");
    yline(0, "HandleVisibility","off", 'Alpha',0.3, 'Color','k', 'LineWidth',3);

    legend("Location", "best", "Orientation","horizontal", "NumColumns",5);

    set(gca,'FontSize', 14);
    xlim(relSegment);
    ylim([-1,1]);
    title(strcat("Planned lane offset - Driver", {' '}, num2str(driverID)));

    subplot(3,1,2);
    plot(X_abs, segment_m(:,indexes.LaneCurvature), 'LineWidth', 2, 'color', 'k');
    grid on; hold on;
    for j=1:length(straightBoundingBoxes)/2
        if (j==1)
            fill([X_abs(straightBoundingBoxes(2*j-1)); X_abs(straightBoundingBoxes(2*j)); X_abs(straightBoundingBoxes(2*j)); X_abs(straightBoundingBoxes(2*j-1))], ...
                [-1;-1;1;1], 'g', 'FaceAlpha', 0.2, 'DisplayName', 'Straight sections', 'LineStyle','none');
        else
            fill([X_abs(straightBoundingBoxes(2*j-1)); X_abs(straightBoundingBoxes(2*j)); X_abs(straightBoundingBoxes(2*j)); X_abs(straightBoundingBoxes(2*j-1))], ...
                [min(segment_m(:,indexes.LaneCurvature)*2);min(segment_m(:,indexes.LaneCurvature)*2);max(segment_m(:,indexes.LaneCurvature)*2);max(segment_m(:,indexes.LaneCurvature)*2)], 'g', 'FaceAlpha', 0.2, 'LineStyle','none', "HandleVisibility","off");
        end
    end
    xlabel('X-UTM(m)'); ylabel("$\kappa(1/m)$");        
    set(gca,'FontSize', 14);
    xlim(relSegment);
    title("Road curvature");
    ylim([-3e-3, 3e-3]);

    subplot(3,1,3);
    leftEdgeCoordinates = [X_abs-segment_m(:,indexes.LaneEdgePositionLeft).*sin(segment_m(:,indexes.theta_calc)) ...
        segment_m(:,indexes.Y_abs)+segment_m(:,indexes.LaneEdgePositionLeft).*cos(segment_m(:,indexes.theta_calc))];
    rightEdgeCoordinates = [X_abs-segment_m(:,indexes.LaneEdgePositionRight).*sin(segment_m(:,indexes.theta_calc)) ...
        segment_m(:,indexes.Y_abs)+segment_m(:,indexes.LaneEdgePositionRight).*cos(segment_m(:,indexes.theta_calc))];
    midLaneCoordinates = [X_abs-segment_m(:,indexes.c0).*sin(segment_m(:,indexes.theta_calc)) ...
        segment_m(:,indexes.Y_abs)+segment_m(:,indexes.c0).*cos(segment_m(:,indexes.theta_calc))];
    fill([leftEdgeCoordinates(:,1); flip(rightEdgeCoordinates(:,1))], [leftEdgeCoordinates(:,2); flip(rightEdgeCoordinates(:,2))], 'g');
    hold on;
    xlim(relSegment);
    for j=1:length(straightBoundingBoxes)/2
        if (j==1)
            fill([X_abs(straightBoundingBoxes(2*j-1)); X_abs(straightBoundingBoxes(2*j)); X_abs(straightBoundingBoxes(2*j)); X_abs(straightBoundingBoxes(2*j-1))], ...
                [min(get(gca,'ylim'));min(get(gca,'ylim'));max(get(gca,'ylim'));max(get(gca,'ylim'))], 'g', 'FaceAlpha', 0.2, 'DisplayName', 'Straight sections', 'LineStyle','none');
        else
            fill([X_abs(straightBoundingBoxes(2*j-1)); X_abs(straightBoundingBoxes(2*j)); X_abs(straightBoundingBoxes(2*j)); X_abs(straightBoundingBoxes(2*j-1))], ...
                [min(get(gca,'ylim'));min(get(gca,'ylim'));max(get(gca,'ylim'));max(get(gca,'ylim'))], 'g', 'FaceAlpha', 0.2, 'LineStyle','none', "HandleVisibility","off");
        end
    end
    plot(X_abs, segment_m(:,indexes.Y_abs), 'color', 'k', 'LineWidth', 2, 'DisplayName', 'Reference');
    plot(midLaneCoordinates(:,1), midLaneCoordinates(:,2), 'color', 'k', 'LineWidth', 2, 'LineStyle', '--', 'DisplayName', 'Centerline');
    grid on;
    xlabel("X-UTM(m)"); ylabel("Y-UTM(m)");
    
    dir = mean(diff(segment_m(:,indexes.X_abs)));
    if (dir < 0)
        annotation('textarrow',[0.7 0.65],[0.22 0.22], 'String','Driving direction', 'FontSize', 14);
    else
        annotation('textarrow',[0.3 0.35],[0.22 0.22], 'String','Driving direction', 'FontSize', 14);
    end
    ylims = get(gca,'ylim');
    set(gca,'FontSize', 14);
    title("Vehicle path");
    ylim(ylims);
end

