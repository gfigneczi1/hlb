function [driver_ID, style] = lanechange_measurement_style(driver_ID, driver_IDs, style, filename, size)
%LANECHANGE_MEASUREMENT_STYLE Summary of this function goes here
%   Detailed explanation goes here
name = strsplit(filename, '.'); name = name{1};
name = strsplit(name, '_'); name = name{end};

if contains(filename, 'style')
    if name == "comfort"
        lc_style = 1;
    elseif name == "dynamic"
        lc_style = 2;
    end
    for i=1:size
        driver_ID{end+1, 1} = 0;
        style{end+1, 1} = lc_style;
    end
else
    ID = getfield(driver_IDs, name);
    for i=1:size
        driver_ID{end+1, 1} = ID;
        style{end+1, 1} = 0;
    end
end
end

