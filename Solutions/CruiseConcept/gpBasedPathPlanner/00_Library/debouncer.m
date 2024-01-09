function [rawData] = debouncer(rawData, dataField, condition, time)
%DEBOUNCER Summary of this function goes here
%   Detailed explanation goes here
dt = 0.01;
time = time / dt; % convert sec to msec
if (isstruct(rawData))
    data_index = find(rawData.(dataField) ~= condition);
    first_ind_mem = 1;
    debounce_index = [];
    for i=1:length(data_index)-1
        if (data_index(i+1) - data_index(i)) ~= 1 ||...
                data_index(i+1) == data_index(end)
            if (i - first_ind_mem < time) &&...
                    data_index(i+1) == data_index(end)
                if (size(data_index,1) > size(data_index,2))
                    debounce_index = [debounce_index data_index(first_ind_mem:i+1)'];
                else
                    debounce_index = [debounce_index data_index(first_ind_mem:i+1)];
                end
            elseif (i - first_ind_mem < time)
                if (size(data_index,1) > size(data_index,2))
                    debounce_index = [debounce_index data_index(first_ind_mem:i)'];
                else
                    debounce_index = [debounce_index data_index(first_ind_mem:i)];
                end
            end
            first_ind_mem = i + 1;
        end
    end
    rawData.(dataField)(debounce_index) = condition;
else
    data_index = find(rawData ~= condition);
    first_ind_mem = 1;
    debounce_index = [];
    for i=1:length(data_index)-1
        if (data_index(i+1) - data_index(i)) ~= 1 ||...
                data_index(i+1) == data_index(end)
            if (i - first_ind_mem < time) &&...
                    data_index(i+1) == data_index(end)
                debounce_index = [debounce_index data_index(first_ind_mem:i+1)];
            elseif (i - first_ind_mem < time)
                debounce_index = [debounce_index data_index(first_ind_mem:i)];
            end
            first_ind_mem = i + 1;
        end
    end
    rawData(debounce_index) = condition;
end

