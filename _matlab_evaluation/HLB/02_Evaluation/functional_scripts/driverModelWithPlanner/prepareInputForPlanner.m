function [segment, segment_m, indexes] = prepareInputForPlanner(segment)

    c0_filtered = movmean(0.5*(segment.LaneEdgePositionLeft+segment.LaneEdgePositionRight), 2);

    if(~isfield(segment, "X_abs"))
          segment.X_abs = segment.LongPos_abs * 40075000 .* cos(segment.LatPos_abs*pi()/180) / 360;
          segment.Y_abs = segment.LatPos_abs * 111.32*1000;
    end
    if (~isfield(segment, "theta_calc"))
        dir = (mean(diff(segment.X_abs))); % +: "north", -: "south"
            if (dir > 0)
                segment.theta_calc = theta_recalc(segment, 1);
            else
                segment.theta_calc = theta_recalc(segment, 0);
            end
    end

    for j=1:size(segment.LaneCurvature,1)
        corridor(j,1:2) = pos_tf2GPS(segment.X_abs(j),segment.Y_abs(j),segment.theta_calc(j),c0_filtered(j));
        map(j,1:2) = pos_tf2GPS(segment.X_abs(j),segment.Y_abs(j),segment.theta_calc(j),segment.LaneEdgePositionLeft(j));
        map(j,3:4) = pos_tf2GPS(segment.X_abs(j),segment.Y_abs(j),segment.theta_calc(j),segment.LaneEdgePositionRight(j));
    end

    segment.corrX = corridor(:,1);
    segment.corrY = corridor(:,2);
    segment.corrLeftX = map(:,1);
    segment.corrLeftY = map(:,2);
    segment.corrRightX = map(:,3);
    segment.corrRightY = map(:,4);
    segment.c0 = 0.5*(segment.LaneEdgePositionLeft+segment.LaneEdgePositionRight);
    % recalculating orientation and curvature
    segment.orientationVector = (pt1_filter(segment.LaneOrientation,1)) + segment.theta_calc;
    segment.curvatureVector = movmean(segment.LaneCurvature,100);
    dT = mean(diff(segment.Relative_time));
    dx = segment.VelocityX*dT;
    
    %% calculating road curvature gradient
    c3 = diff(segment.LaneCurvature*2)./dx(2:end);
    c3 = [c3; c3(end)];
    c3 = movmean(c3,200);
    
    segment.c3 = c3;
    clear corridor;

    % SAVING out field indices
    varnames = segment.Properties.VariableNames;
    [~, indexes.X_abs] = ismember('X_abs', varnames);
    [~, indexes.Y_abs] = ismember('Y_abs', varnames);
    [~, indexes.theta_calc] = ismember('theta_calc', varnames);
    [~, indexes.c0] = ismember('c0',varnames);
    [~, indexes.LaneEdgePositionLeft] = ismember('LaneEdgePositionLeft',varnames);
    [~, indexes.LaneEdgePositionRight] = ismember('LaneEdgePositionRight',varnames);
    [~, indexes.LaneOrientation] = ismember('LaneOrientation',varnames);
    [~, indexes.LaneCurvature] = ismember('LaneCurvature',varnames);
    [~, indexes.c3] = ismember('c3',varnames);
    [~, indexes.orientationVector] = ismember('orientationVector',varnames);
    [~, indexes.curvatureVector] = ismember('curvatureVector',varnames);
    [~, indexes.corrX] = ismember('corrX',varnames);
    [~, indexes.corrY] = ismember('corrY',varnames);
    [~, indexes.corrLeftX] = ismember('corrLeftX',varnames);
    [~, indexes.corrLeftY] = ismember('corrLeftY',varnames);
    [~, indexes.corrRightX] = ismember('corrRightX',varnames);
    [~, indexes.corrRightY] = ismember('corrRightY',varnames);
    [~, indexes.VelocityX] = ismember('VelocityX',varnames);
    [~, indexes.AccelerationX] = ismember('AccelerationX', varnames);
    [~, indexes.Acceleration_Y] = ismember('Acceleration_Y', varnames);
    [~, indexes.Relative_time] = ismember('Relative_time',varnames);
    [~, indexes.YawRate] = ismember('YawRate',varnames);
    [~, indexes.OncomingTrafficType] = ismember('OncomingTrafficType', varnames);
    [~, indexes.FrontTrafficType] = ismember('FrontTrafficType', varnames);
    [~, indexes.OncomingVehicleTimeToPass] = ismember('OncomingVehicleTimeToPass', varnames);
    [~, indexes.drID] = ismember('drID', varnames);

    segment_m = table2array(segment);

    if (indexes.OncomingTrafficType == 0)
        segment_m(:, end+1) = 0;
        indexes.OncomingTrafficType = size(segment_m,2);
    end
   
    % adding extra channels
    segment_m(:,end+1) = movmean(segment_m(:, indexes.YawRate) ./ segment_m(:, indexes.VelocityX), 25);
    indexes.vehicleCurvature = size(segment_m,2);
    segment_m(:,end) = segment_m(:,end) - mean(segment_m(:,end));
    segment_m(1:end-1,end+1) = movmean(diff(segment_m(:,indexes.vehicleCurvature))./diff(segment_m(:, indexes.Relative_time)),25);
    segment_m(end,end) = segment_m(end-1,end);
    segment_m(:,end) = segment_m(:,end) - mean(segment_m(:,end));
    indexes.vehicleCurvatureChange = size(segment_m,2);

    segment_m(:,end+1) =  movmean(segment_m(:, indexes.YawRate) .* segment_m(:, indexes.VelocityX) - (segment_m(:, indexes.VelocityX).^2) .* segment_m(:, indexes.LaneCurvature)*2, 25);
    indexes.ayRel = size(segment_m,2);
    segment_m(:,end) = segment_m(:,end) - mean(segment_m(:,end));
    segment_m(1:end-1,end+1) = movmean(diff(segment_m(:,indexes.ayRel))./diff(segment_m(:, indexes.Relative_time)), 25);
    segment_m(end,end) = segment_m(end-1,end);
    segment_m(:,end) = segment_m(:,end) - mean(segment_m(:,end));
    indexes.jyRel = size(segment_m,2);
    
    % predicted time to lane crossing (TTLC)
    segment_m(1:end-1,end+1) = diff(-segment_m(:,indexes.c0))./diff(segment_m(:,indexes.Relative_time));
    segment_m(end,end) = segment_m(end-1,end);
    segment_m(:,end) = movmean(segment_m(:,end),10);
    indexes.vy = size(segment_m,2);
    segment_m(:,end+1) = segment_m(:,indexes.VelocityX).*sin((-atan(segment_m(:,indexes.LaneOrientation)))) + segment_m(:,indexes.vy).*cos((-atan(segment_m(:,indexes.LaneOrientation))));
    indexes.veta = size(segment_m,2);
    segment_m(:,end+1) = segment_m(:,indexes.AccelerationX).*sin((-atan(segment_m(:,indexes.LaneOrientation)))) + segment_m(:,indexes.Acceleration_Y).*cos((-atan(segment_m(:,indexes.LaneOrientation))));
    indexes.aeta = size(segment_m,2);

    driftingSide = segment_m(:,indexes.veta)>0; % 1: left side is critical, 0: right side is critical
    dcrit = zeros(size(segment_m,1),1);
    dcrit(driftingSide) = segment_m(driftingSide,indexes.LaneEdgePositionLeft);
    dcrit(~driftingSide) = segment_m(~driftingSide,indexes.LaneEdgePositionRight);
    ttcl1 = (-segment_m(:,indexes.veta)+(segment_m(:,indexes.veta).^2+2*segment_m(:,indexes.aeta).*dcrit).^0.5)./segment_m(:,indexes.aeta);
    ttcl2 = (-segment_m(:,indexes.veta)-(segment_m(:,indexes.veta).^2+2*segment_m(:,indexes.aeta).*dcrit).^0.5)./segment_m(:,indexes.aeta);
    % selecting the TTCL-s from the vectors
    relevantttcl = (imag(ttcl1)==0 & ttcl1>=0) | (imag(ttcl2)==0 & ttcl2>=0);
    ttcl1(imag(ttcl1)~=0 | real(ttcl1)<0) = inf;
    ttcl2(imag(ttcl2)~=0 | real(ttcl2)<0) = inf;

    TTCL = min(ttcl1, ttcl2);
    TTCL(~relevantttcl) = inf;

    segment_m(:,end+1) = TTCL;
    indexes.TTCL = size(segment_m,2);

    % average path (needed for lane wandering analysis)
    dT = mean(diff(segment_m(:,indexes.Relative_time)));
    window = floor(9/dT);
    c0_filtered = movmean(segment_m(:,indexes.c0), window);
    
    for j=1:size(segment.LaneCurvature,1)
        modifiedReference(j,1:2) = pos_tf2GPS(segment_m(j,indexes.corrX),segment_m(j,indexes.corrY),segment_m(j, indexes.orientationVector),-c0_filtered(j));
    end

    segment_m(:,end+1) = modifiedReference(:,1);
    indexes.X_abs_mod = size(segment_m,2);
    segment_m(:,end+1) = modifiedReference(:,2);
    indexes.Y_abs_mod = size(segment_m,2);
    segment_m(:,end+1) = c0_filtered;
    indexes.c0_filtered = size(segment_m,2);

end

