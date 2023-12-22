function [slowOffset, offsetError, laneOffset] = functional_calculateStraightLineOffset(segment_m, indexes, fc, dT)
    % butterworth filter is applied
    fs = 1/dT;

    [b,a] = butter(2,fc/(fs/2));
    laneOffset = -segment_m(:, indexes.c0); % negative sign due to different coordinate system
    slowOffset = filter(b,a,-segment_m(:, indexes.c0));
    offsetError = -segment_m(:,indexes.c0)-slowOffset;
end

