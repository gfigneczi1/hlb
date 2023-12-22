function [driver_ID] = create_ID(matFiles)
%CREATE_ID Summary of this function goes here
%   Detailed explanation goes here
names = double.empty;
for i=1:length(matFiles)
    filename = matFiles(i).name;
    if ~contains(filename, 'style')
        name = strsplit(filename, '.'); name = name{1};
        name = strsplit(name, '_'); name = name{end};
        names = [names string(name)];
    end
end
names = unique(names);
if ~isempty(names)
    for i=1:length(names)
        driver_ID.(names(i)) = i;
    end
else
    driver_ID = double.empty;
end
end

