function [orientation, curvature] = calcPathGeometry(path)
    % Orientation: the orientation angle function in radians
    % Curvature: the smoothed curvature of the path in 1/m
    % derivation of driving curvature
    orientation = [diff(path(:,2))./diff(path(:,1)); 0];
    curvature = [diff(orientation)./diff(path(:,1)); 0];
    curvature = movmean(curvature,50);
end

