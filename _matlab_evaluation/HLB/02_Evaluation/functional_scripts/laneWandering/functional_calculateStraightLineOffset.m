function [slowOffset, offsetError, laneOffset] = functional_calculateStraightLineOffset(segment_m, indexes, fc, dT)
    % butterworth filter is applied
    fs = 1/dT;
    
    window = round(fs/fc);

    %[b,a] = butter(2,fc/(fs/2));
    %[b,a] = cheby1(2,10,fc/(fs/2));
    laneOffset = -segment_m(:, indexes.c0); % negative sign due to different coordinate system
    slowOffset = movmean(-segment_m(:, indexes.c0), 180);
    %slowOffset = filter(b,a,-segment_m(:, indexes.c0));
    offsetError = -segment_m(:,indexes.c0)-slowOffset;    
end

