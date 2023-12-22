function [rawDatas, lapnumber] = circle_cutter(rawData)
X_abs = rawData.LongPos_abs' * 40075000 .* cos(rawData.LatPos_abs'*pi()/180) / 360;
Y_abs = rawData.LatPos_abs' * 111.32*1000;

X_start_Rural_south = [X_abs(1)-100, X_abs(1)+100, X_abs(1)-100, X_abs(1)+100];
Y_start_Rural_south = [Y_abs(1)-100, Y_abs(1)-100, Y_abs(1)+100, Y_abs(1)+100];

laps = (inpolygon(X_abs, Y_abs, X_start_Rural_south,...
        Y_start_Rural_south));

rawDatas = {};

lapstarts = find(diff(laps)==1);
if (laps(1) == 1)
    lapstarts = [1; lapstarts];
end

if (numel(lapstarts) == 1)
    rawDatas{1} = rawData;
else
    rawData.lapnumber = zeros(size(X_abs,1),1)';
    j = 1;
    for i=1:numel(lapstarts)-1
        rawData.lapnumber(lapstarts(i):lapstarts(i+1)-1) = i;
        if (lapstarts(i+1)-lapstarts(i)) > 500
            rawDatas{j} = structure_filter(rawData,"lapnumber",i);
            j=j+1;
        end
    end
    if (numel(lapstarts) == 1 || isempty(rawDatas))
        rawDatas = {};
    end
end
lapnumber = numel(rawDatas);

end

