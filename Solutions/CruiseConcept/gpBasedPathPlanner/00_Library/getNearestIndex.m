function [nearestIndex] = getNearestIndex(X,Y,P)
    distances = ((X-P(1)).^2+(Y-P(2)).^2).^0.5;
    nearestIndex = find(distances==min(distances),1);
end

