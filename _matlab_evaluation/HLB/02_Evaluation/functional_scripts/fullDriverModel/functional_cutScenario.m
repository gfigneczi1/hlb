function scenario = functional_cutScenario(segment_m, indexes, nearestIndex, vehicleState, window)
    % window = [pastDistance FutureDistance]
    vehiclePose = [vehicleState.X vehicleState.Y vehicleState.theta];
    distances = ((segment_m(:, indexes.X_abs) - vehiclePose(1)).^2+(segment_m(:, indexes.Y_abs) - vehiclePose(2)).^2).^0.5;
    if (window(1) == 0)
        scenarioStart = nearestIndex;
    else
        scenarioStart = find(distances(1:nearestIndex) <= window(1),1);
    end
    scenarioStop = find(distances(nearestIndex:end) > window(2),1);
    if (isempty(scenarioStart) || isempty(scenarioStop))
        scenario = [];
    else
        scenarioStop = nearestIndex + scenarioStop - 1;
        if (any(abs(diff(distances(scenarioStart:scenarioStop))) > 3))
            % there is discontinuity in the data
            scenario = [];
        else
            scenario = segment_m(scenarioStart:scenarioStop, :);
        end
    end
end


