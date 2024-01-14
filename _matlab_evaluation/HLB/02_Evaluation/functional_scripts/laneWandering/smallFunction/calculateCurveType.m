function curveTypes = calculateCurveType(segment_m, indexes, name)

global temp_folder_path plots_folder_name
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

    g = figure();

    plot(segment_m(:,indexes.X_abs), segment_m(:, indexes.Y_abs));
    grid on; hold on;
    plot(segment_m(curveTypes==1, indexes.X_abs), segment_m(curveTypes==1, indexes.Y_abs), 'ko');
    legend('original', 'straights');
    savefig(g, fullfile(temp_folder_path, plots_folder_name,...
                strcat('Map', name(1:end-4), '.fig')));
    saveas(g, fullfile(temp_folder_path, plots_folder_name,...
                strcat('Map_', name(1:end-4), '.png')));
    close(g);

end

