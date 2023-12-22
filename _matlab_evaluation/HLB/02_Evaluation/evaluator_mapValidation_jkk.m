function evaluator_mapValidation_jkk (segments,config)
    % Checking input data validity
    lookAheadDistance = 100; % meters
    for i=1:size(segments.segments,2)
        
        segment = segments.segments(i).segment;
        segment.X_abs = segment.LongPos_abs * 40075000 .* cos(segment.LatPos_abs*pi()/180) / 360;
        segment.Y_abs = segment.LatPos_abs * 111.32*1000;

        segment.Left_Index = zeros(size(segment.X_abs,1));
        segment.Right_Index = zeros(size(segment.X_abs,1));
        mapLeft = segments.maps.(segments.segments(i).lane)(:,1:2);
        mapRight = segments.maps.(segments.segments(i).lane)(:,3:4);
        segment.theta_calc = movmean(segment.theta_calc,100);
        % Looping through all vehicle positions and transforming it to the
        % ego frame
        for j = 1:size(segment,1)
            mapLeftLocal = GPS2ego([segment.X_abs(j),segment.Y_abs(j), segment.theta_calc(j)],mapLeft);
            mapRightLocal = GPS2ego([segment.X_abs(j),segment.Y_abs(j), segment.theta_calc(j)],mapRight);
            % finding the right subsegment of the maps
            distances = (mapLeftLocal(:,1).^2 + mapLeftLocal(:,2).^2).^0.5;
            pointsWithinLookAheadDistanceLeft = mapLeftLocal(distances < lookAheadDistance,:);
            distances = (mapRightLocal(:,1).^2 + mapRightLocal(:,2).^2).^0.5;
            pointsWithinLookAheadDistanceRight = mapRightLocal(distances < lookAheadDistance,:);
            if (size(pointsWithinLookAheadDistanceLeft,1) < 2 || size(pointsWithinLookAheadDistanceRight,1) < 2)
                segment.c01_left(j) = 0;
                segment.c01_right(j) = 0;
                segment.c1(j) = 0;
                segment.c2(j) = 0;
            else
                % smoothing before polyline fit - only when closest point
                % is close enough
                if (pointsWithinLookAheadDistanceLeft(1,1) < 5 && pointsWithinLookAheadDistanceRight(1,1) < 5)
                    pointsWithinLookAheadDistanceLeft = [linspace(pointsWithinLookAheadDistanceLeft(1,1),pointsWithinLookAheadDistanceLeft(end,1),500)' spline(pointsWithinLookAheadDistanceLeft(:,1), pointsWithinLookAheadDistanceLeft(:,2), linspace(pointsWithinLookAheadDistanceLeft(1,1),pointsWithinLookAheadDistanceLeft(end,1),500))']; 
                    pointsWithinLookAheadDistanceRight = [linspace(pointsWithinLookAheadDistanceRight(1,1),pointsWithinLookAheadDistanceRight(end,1),500)' spline(pointsWithinLookAheadDistanceRight(:,1), pointsWithinLookAheadDistanceRight(:,2), linspace(pointsWithinLookAheadDistanceRight(1,1),pointsWithinLookAheadDistanceRight(end,1),500))']; 
                    cleft = polyfit(pointsWithinLookAheadDistanceLeft(:,1),pointsWithinLookAheadDistanceLeft(:,2),2);
                    cright = polyfit(pointsWithinLookAheadDistanceRight(:,1),pointsWithinLookAheadDistanceRight(:,2),2);
                    segment.c01_left(j) = cleft(3);
                    segment.c01_right(j) = cright(3);
                    segment.c1(j) = atan((cleft(2)+cright(2)/2));
                    segment.c2(j) = cleft(1)+cright(1);
                else
                    segment.c01_left(j) = 0;
                    segment.c01_right(j) = 0;
                    segment.c1(j) = 0;
                    segment.c2(j) = 0;
                end
            end
        end
        segment = table2struct(segment,"ToScalar",true);
        signals = fieldnames(segment);
        for j=1:length(signals)
            fieldname = convertCharsToStrings(signals(j));
            segment.(fieldname) = segment.(fieldname)';
        end
        save (fullfile(config.root,strcat(segments.segments(i).name,"_extended.mat")),'segment');

    end