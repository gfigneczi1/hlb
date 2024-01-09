function [subsegment] = cutSubsegment(segment, startIndex, endIndex)
datafield = fieldnames(segment);
datafield = datafield(1:end-3);

for i = 1:numel(datafield)
    rawData_filtered.(datafield{i}) = segment.(datafield{i})(startIndex:endIndex);
end
subsegment = rawData_filtered;
end

