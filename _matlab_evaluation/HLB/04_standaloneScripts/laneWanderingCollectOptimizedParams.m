close all;
clear all;

paramFilesLeft = dir(fullfile("C:\Users\igg2bp\Documents\00_Projects\Research\laneWandering\_dataIn\optimizedParameters", "\*.mat"));
drID_prev = "";
for i=1:length(paramFilesLeft)
    drID = convertCharsToStrings(paramFilesLeft(i).name(1:find(paramFilesLeft(i).name=='.',1)-1));
    if (drID~=drID_prev && drID_prev~="")
        save(fullfile(paramFilesLeft(i).folder, strcat(name_prev(1:find(paramFilesLeft(i).name=='.',1)-1), '.mat')), 'OptimizedParameters');
        OptimizedParameters = [];
    end
    snippetID = str2num(paramFilesLeft(i).name(length(paramFilesLeft(i).name)-find(paramFilesLeft(i).name(end:-1:1)=='_',1)+2:end-4));

    p = load(fullfile(paramFilesLeft(i).folder, paramFilesLeft(i).name));
    OptimizedParameters(snippetID,:) = p.y;
    
    drID_prev = drID;
    name_prev = paramFilesLeft(i).name;
end