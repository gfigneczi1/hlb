function [segment, segment_m, indexes] = prepareInputForPlanner(segment)

    c0_filtered = movmean(0.5*(segment.c01_left+segment.c01_right), 2);


    for j=1:size(segment.c2,1)
        corridor(j,1:2) = pos_tf2GPS(segment.X_abs(j),segment.Y_abs(j),segment.theta_calc(j),c0_filtered(j));
        map(j,1:2) = pos_tf2GPS(segment.X_abs(j),segment.Y_abs(j),segment.theta_calc(j),segment.c01_left(j));
        map(j,3:4) = pos_tf2GPS(segment.X_abs(j),segment.Y_abs(j),segment.theta_calc(j),segment.c01_right(j));
    end

    segment.corrX = corridor(:,1);
    segment.corrY = corridor(:,2);
    segment.corrLeftX = map(:,1);
    segment.corrLeftY = map(:,2);
    segment.corrRightX = map(:,3);
    segment.corrRightY = map(:,4);
    segment.c0 = 0.5*(segment.c01_left+segment.c01_right);
    % recalculating orientation and curvature
    segment.orientationVector = atan(pt1_filter(segment.c1,1)) + segment.theta_calc;
    segment.curvatureVector = movmean(segment.c2,100);
    segment.c3 = zeros(size(segment.c0));
    clear corridor;

    % SAVING out field indices
    varnames = segment.Properties.VariableNames;
    [~, indexes.X_abs] = ismember('X_abs', varnames);
    [~, indexes.Y_abs] = ismember('Y_abs', varnames);
    [~, indexes.theta_calc] = ismember('theta_calc', varnames);
    [~, indexes.c0] = ismember('c0',varnames);
    [~, indexes.c01_left] = ismember('c01_left',varnames);
    [~, indexes.c01_right] = ismember('c01_right',varnames);
    [~, indexes.c1] = ismember('c1',varnames);
    [~, indexes.c2] = ismember('c2',varnames);
    [~, indexes.c3] = ismember('c3',varnames);
    [~, indexes.orientationVector] = ismember('orientationVector',varnames);
    [~, indexes.curvatureVector] = ismember('curvatureVector',varnames);
    [~, indexes.corrX] = ismember('corrX',varnames);
    [~, indexes.corrY] = ismember('corrY',varnames);
    [~, indexes.corrLeftX] = ismember('corrLeftX',varnames);
    [~, indexes.corrLeftY] = ismember('corrLeftY',varnames);
    [~, indexes.corrRightX] = ismember('corrRightX',varnames);
    [~, indexes.corrRightY] = ismember('corrRightY',varnames);
    [~, indexes.velocityX] = ismember('VelocityX_ESP',varnames);
    [~, indexes.accelerationX] = ismember('AccelerationX_ESP', varnames);
    [~, indexes.accelerationY] = ismember('AccelerationY_ESP', varnames);
    [~, indexes.q_T0] = ismember('q_T0',varnames);
    [~, indexes.yawRate] = ismember('yawRateESP',varnames);
    [~, indexes.oncomingTraffic] = ismember('oncomingTraffic', varnames);
    [~, indexes.drID] = ismember('drID', varnames);

    segment_m = table2array(segment);

    if (indexes.oncomingTraffic == 0)
        segment_m(:, end+1) = 0;
        indexes.oncomingTraffic = size(segment_m,2);
    end
   
    % adding extra channels
    segment_m(:,end+1) = movmean(segment_m(:, indexes.yawRate) ./ segment_m(:, indexes.velocityX), 25);
    indexes.vehicleCurvature = size(segment_m,2);
    segment_m(:,end) = segment_m(:,end) - mean(segment_m(:,end));
    segment_m(1:end-1,end+1) = movmean(diff(segment_m(:,indexes.vehicleCurvature))./diff(segment_m(:, indexes.q_T0)),25);
    segment_m(end,end) = segment_m(end-1,end);
    segment_m(:,end) = segment_m(:,end) - mean(segment_m(:,end));
    indexes.vehicleCurvatureChange = size(segment_m,2);

    segment_m(:,end+1) =  movmean(segment_m(:, indexes.yawRate) .* segment_m(:, indexes.velocityX) - (segment_m(:, indexes.velocityX).^2) .* segment_m(:, indexes.c2)*2, 25);
    indexes.ayRel = size(segment_m,2);
    segment_m(:,end) = segment_m(:,end) - mean(segment_m(:,end));
    segment_m(1:end-1,end+1) = movmean(diff(segment_m(:,indexes.ayRel))./diff(segment_m(:, indexes.q_T0)), 25);
    segment_m(end,end) = segment_m(end-1,end);
    segment_m(:,end) = segment_m(:,end) - mean(segment_m(:,end));
    indexes.jyRel = size(segment_m,2);
    
    % predicted time to lane crossing (TTLC)
    segment_m(1:end-1,end+1) = diff(-segment_m(:,indexes.c0))./diff(segment_m(:,indexes.q_T0));
    segment_m(end,end) = segment_m(end-1,end);
    segment_m(:,end) = movmean(segment_m(:,end),10);
    indexes.vy = size(segment_m,2);
    segment_m(:,end+1) = segment_m(:,indexes.velocityX).*sin((-atan(segment_m(:,indexes.c1)))) + segment_m(:,indexes.vy).*cos((-atan(segment_m(:,indexes.c1))));
    indexes.veta = size(segment_m,2);
    segment_m(:,end+1) = segment_m(:,indexes.accelerationX).*sin((-atan(segment_m(:,indexes.c1)))) + segment_m(:,indexes.accelerationY).*cos((-atan(segment_m(:,indexes.c1))));
    indexes.aeta = size(segment_m,2);

    driftingSide = segment_m(:,indexes.veta)>0; % 1: left side is critical, 0: right side is critical
    dcrit = zeros(size(segment_m,1),1);
    dcrit(driftingSide) = segment_m(driftingSide,indexes.c01_left);
    dcrit(~driftingSide) = segment_m(~driftingSide,indexes.c01_right);
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
    dT = mean(diff(segment_m(:,indexes.q_T0)));
    window = floor(9/dT);
    c0_filtered = movmean(segment_m(:,indexes.c0), window);
    
    for j=1:size(segment.c2,1)
        modifiedReference(j,1:2) = pos_tf2GPS(segment_m(j,indexes.corrX),segment_m(j,indexes.corrY),segment_m(j, indexes.orientationVector),-c0_filtered(j));
    end

    segment_m(:,end+1) = modifiedReference(:,1);
    indexes.X_abs_mod = size(segment_m,2);
    segment_m(:,end+1) = modifiedReference(:,2);
    indexes.Y_abs_mod = size(segment_m,2);

end

