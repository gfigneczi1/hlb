function curveTypes = calculateCurveType(segment_m, indexes)

    % return: curveTypes: 0 = unknown, 1=straight, 2=left, 3=right
    thd = 3.5e-4;
    curveTypes = zeros(size(segment_m,1),1);
    curvature = movmean(segment_m(:,indexes.c2)*2, 50);
    straightLine = abs(curvature) < thd;
    straightLine = morphologyOpen(straightLine, 100);
    straightLine = morphologyOpen(straightLine, 50);
    straightLine = morphologyClose(straightLine, 100);
    curveTypes(straightLine) = 1;
    
    leftCurve = curvature>=thd;
    leftCurve = morphologyClose(leftCurve,50);
    rightCurve = curvature<=-thd;
    rightCurve = morphologyClose(rightCurve,50);
    curveTransition = leftCurve&rightCurve;
    curveTypes(leftCurve&~curveTransition) = 2;
    curveTypes(rightCurve&~curveTransition) = 3;
    curveTypes(curveTransition) = 2.5;
end
