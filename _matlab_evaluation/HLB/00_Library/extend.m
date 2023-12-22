function [rawData] = extend(rawData, dataField, condition, time)
dt = 0.05;
time = time / dt;
if (isstruct(rawData))
    data_index = find(diff(rawData.(dataField)) >0);
    if (time < 0)
        time = -time;
        for i=1:length(data_index)
            rawData.(dataField)(max(1, data_index(i)-round(time)):data_index(i)) = condition;
        end
    else
        for i=1:length(data_index)
            rawData.(dataField)(data_index(i):min(length(rawData.(dataField)),data_index(i)+round(time))) = condition;
        end
    end
else
    
end
end

