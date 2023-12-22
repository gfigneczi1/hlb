clear;
close all;

input_folder = "../../_temp";

% load mat files
matFiles = dir(fullfile(input_folder,'/*.mat'));
labelFiles = dir(fullfile(input_folder,'/*.xlsx'));
for fileID = 1:length(matFiles)
    rawData = load(fullfile(matFiles(fileID).folder,matFiles(fileID).name));
    drID = matFiles(fileID).name(1:5);
    for j=1:length(labelFiles)
        if (labelFiles(j).name(1:5) == drID)
            labels = readtable(fullfile(matFiles(j).folder,labelFiles(j).name));
            break;
        end
    end
    % adding extra channel of oncoming traffic 
    time = rawData.segment.q_T0 - rawData.segment.q_T0(1);
    oncomingTraffic = zeros(1,size(rawData.segment.q_T0,2));
    for i=1:size(labels,1)
        % mapping: 1 - truck, 2 - normal vehicle, 2x - convoy (21 - convoy
        % truck, 22 - convoy normal, 23 - convoy-mixed, 
        % 0 - no oncoming traffic
        t0 = find(time>=labels.trafficStart_dec(i),1);
        t1 = find(time>=labels.trafficEnd_dec(i),1);
        if (contains(labels.objectType(i), "overtake"))
            oncomingTraffic(t0:t1) = 30;
        elseif (contains(labels.objectType(i), "convoy") && contains(labels.objectType(i), "truck"))
            oncomingTraffic(t0:t1) = 21;
        elseif (contains(labels.objectType(i), "convoy") && contains(labels.objectType(i), "normal"))
            oncomingTraffic(t0:t1) = 22;
        elseif (contains(labels.objectType(i), "convoy") && contains(labels.objectType(i), "mixed"))
            oncomingTraffic(t0:t1) = 23;
        elseif (contains(labels.objectType(i), "truck"))
            oncomingTraffic(t0:t1) = 1;
        elseif (contains(labels.objectType(i), "normal"))
            oncomingTraffic(t0:t1) = 2;
        else
            oncomingTraffic(t0:t1) = 2;
        end
    end
        
    rawData.segment.oncomingTraffic = oncomingTraffic;

    save(fullfile(matFiles(fileID).folder,matFiles(fileID).name),'-struct', "rawData");
end
