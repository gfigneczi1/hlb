function [rawData] = monotonyCheck(rawData, fieldName)
% Finding extremum
monotony = diff(rawData.(fieldName)) > 0;
monotony(monotony==0) = -1;
rawData.monotony = monotony;
end

