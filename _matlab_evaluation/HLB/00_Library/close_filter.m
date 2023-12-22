function [rawData] = close_filter(rawData)
if (rawData(1) == 1)
    bridges = find(diff(rawData)==1) - find(diff(rawData)==-1);
    r = find(diff(rawData)==-1);
elseif (rawData(end) == 1)
    r = [1; find(diff(rawData)==-1)];
    bridges = find(diff(rawData)==1) - r;
else
    r = find(diff(rawData)==1);
    bridges = find(diff(rawData)==1) - find(diff(rawData)==-1);
end

for i=1:numel(bridges)
    if bridges(i) < 2000
        rawData(r(i):r(i)+bridges(i)) = 1;
    end
end
end

