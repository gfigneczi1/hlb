function extrema = calculateExtrema(offsetError,relevantPoints)
    thd = std(offsetError(relevantPoints==1))*0.5;
    interactionPoints = abs(offsetError) > thd;
    interactionPoints(~relevantPoints) = false;
    interactionPoints = morphologyOpen(interactionPoints, 10);
    interactionPoints = morphologyClose(interactionPoints, 10);
    extrema = findExtremumPoints(offsetError, interactionPoints);
end

