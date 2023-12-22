function [rawData] = redundantPointsFilter(rawData, type)
    datafield = fieldnames(rawData);
    rawData_filtered = rawData;
    condition = 0 * (rawData.GPS_time);
    % diff(rawData.GPS_time == 0);
    if (type == "Time")
        condition(2:end) = diff(rawData.GPS_time) == 0;
    elseif (type == "Position")
        condition(2:end) = (diff(rawData.LongPos_abs) == 0) .* (diff(rawData.LatPos_abs) == 0);
    else
        condition(2:end) = 1;
    end
    for i = 1:numel(datafield)
        rawData_filtered.(datafield{i}) = rawData.(datafield{i})( condition == 0 );
    end
    rawData = rawData_filtered;
end

