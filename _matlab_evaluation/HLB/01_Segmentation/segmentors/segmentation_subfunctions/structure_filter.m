function [rawData, cuttingIndexes] = structure_filter(rawData, dataField, condition, type)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
datafield = fieldnames(rawData);
if (istable(rawData))
    datafield = datafield(1:end-3);
else
    rawData_filtered = rawData;
end
try
    cuttingIndexes(:,1) = find(diff(rawData.(dataField) == condition) == 1)+2;
    cuttingIndexes(:,2) = find(diff(rawData.(dataField) == condition) == -1)+2;
catch
    cuttingIndexes = [];
end
for i = 1:numel(datafield)
    if type == "GT"
        rawData_filtered.(datafield{i}) = rawData.(datafield{i})(rawData.(dataField) > condition);
    elseif type == "EQ"
        rawData_filtered.(datafield{i}) = rawData.(datafield{i})(rawData.(dataField) == condition);
    end
end
rawData = rawData_filtered;
end

