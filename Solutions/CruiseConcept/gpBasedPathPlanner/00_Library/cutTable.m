function [structure] = cutTable(table, indeces, preIndex, postIndex)
% This function returns a structure with n fields, n = number of indeces;
% the cut structure elements are the cut sections of entire dataset(table)
% indeces must be an n x 2 array, where column 1 is startIndex, column 2 is
% endindex

if (istable(table) && ~isempty(indeces))
    datafield = fieldnames(table);
    datafield = datafield(1:end-3);
    n = length(table.(datafield{1}));
    for i=1:size(indeces,1)
        for j = 1:numel(datafield)
            table_filtered.(datafield{j}) = table.(datafield{j})(max(1,indeces(i,1)-preIndex):min(indeces(i,2)+postIndex,n));
        end
        structure(i) = table_filtered;
        clear table_filtered;
    end
else
    structure = struct();
end

end

