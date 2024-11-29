clear;
pathToCsv = "C:\git\KDP\HLB_for_AV_MotionControl\02_Results\MassMeasurements\Dr022\Dr022_trafficLabelsSecs.xlsx";
pathToMat = 'C:\git\KDP\HLB_for_AV_MotionControl\02_Results\MassMeasurements\Dr022\Dr022_62A_2024-10-17.mat_extended.mat';

%% Mapping of object types
% 1 - normal, convoy-normal
% 2 - truck(s), convoy trucks
% 3 - convoy, mixed

% start time: when object definitely appears on video. Note: if there is
% traffic in front, this indicates lower distance-to-appear.

% levelOfPresence: signal which indicates how far the given oncoming object/convoy
% is. This is a calculated signal, assuming nominal speed levels+margin and distance of
% appearance. 

%% Calculation

traffic = readtable(pathToCsv);
load(pathToMat);

oncomingTraffic = zeros(1,size(segment.q_T0,2));
trafficInFront = zeros(1,size(segment.q_T0,2));
timeToPass = zeros(1,size(segment.q_T0,2));

segment.q_T0 = segment.q_T0 - segment.q_T0(1);

for i=1:size(traffic,1)
    if(contains(convertCharsToStrings(traffic.objType{i}), "normal") | isempty(traffic.objType{i}))
        startTime = find(segment.q_T0 > traffic.trafficStart(i), 1);
        endTime = find(segment.q_T0 > traffic.trafficEnd(i), 1);
        oncomingTraffic(startTime:endTime) = 1;
    elseif(contains(convertCharsToStrings(traffic.objType{i}), "convoy"))
        startTime = find(segment.q_T0 > traffic.trafficStart(i), 1);
        endTime = find(segment.q_T0 > traffic.trafficEnd(i), 1);
        oncomingTraffic(startTime:endTime) = 3;
    elseif(contains(convertCharsToStrings(traffic.objType{i}), "truck"))
        startTime = find(segment.q_T0 > traffic.trafficStart(i), 1);
        endTime = find(segment.q_T0 > traffic.trafficEnd(i), 1);
        oncomingTraffic(startTime:endTime) = 2;
    elseif(contains(convertCharsToStrings(traffic.objType{i}), "in overtake"))
        startTime = find(segment.q_T0 > traffic.trafficStart(i), 1);
        endTime = find(segment.q_T0 > traffic.trafficEnd(i), 1);
        oncomingTraffic(startTime:endTime) = 4;
    end
    if(contains(convertCharsToStrings(traffic.objInFrontType{i}), "normal"))
        startTime = find(segment.q_T0 > traffic.trafficInFrontStart(i), 1);
        endTime = find(segment.q_T0 > traffic.trafficInFrontEnd(i), 1);
        trafficInFront(startTime:endTime) = 1;
    elseif(contains(convertCharsToStrings(traffic.objInFrontType{i}), "truck"))
        startTime = find(segment.q_T0 > traffic.trafficInFrontStart(i), 1);
        endTime = find(segment.q_T0 > traffic.trafficInFrontEnd(i), 1);
        trafficInFront(startTime:endTime) = 2;
    end
end
i = 2;
TTA_norm = 2.38;
TTA_truck = 3.38;
while(i<=length(oncomingTraffic))
    if (oncomingTraffic(i) == 1 && oncomingTraffic(i-1) == 0)
        if (isempty(find(oncomingTraffic(i:end)==0,1)))
            i = length(oncomingTraffic)+1;
        else
            timeOfAppearance = segment.q_T0(i);
            if(trafficInFront(i)==0)
                timeOfPass = timeOfAppearance + TTA_norm;
            elseif(trafficInFront(i)==1)
                timeOfPass = timeOfAppearance + TTA_norm*2/3;
            else
                timeOfPass = timeOfAppearance + TTA_norm/2;
            end
            timeToPass(i:i+find(segment.q_T0(i:end) >= timeOfPass,1)-1) = linspace(timeOfPass-timeOfAppearance, 0, find(segment.q_T0(i:end) >= timeOfPass,1));
        end
        i = i+find(oncomingTraffic(i:end)==0,1);
    elseif (oncomingTraffic(i) >= 2 && oncomingTraffic(i-1) == 0)
        if (isempty(find(oncomingTraffic(i:end)==0,1)))
            i = length(oncomingTraffic)+1;
        else
            timeOfAppearance = segment.q_T0(i);
            if(trafficInFront(i)==0)
                timeOfPass = timeOfAppearance + TTA_truck;
            elseif(trafficInFront(i)==1)
                timeOfPass = timeOfAppearance + TTA_truck*2/3;
            else
                timeOfPass = timeOfAppearance + TTA_truck/2;
            end
            timeToPass(i:i+find(segment.q_T0(i:end) >= timeOfPass,1)-1) = linspace(timeOfPass-timeOfAppearance, 0, find(segment.q_T0(i:end) >= timeOfPass,1));
        end
        i = i+find(oncomingTraffic(i:end)==0,1);
    else
        if(trafficInFront(i)==0)
            timeToPass(i) = TTA_norm;
        elseif(trafficInFront(i)==1)
            timeToPass(i) = TTA_norm*2/3;
        else
            timeToPass(i) = TTA_norm/2;
        end
        i=i+1;
    end
end

segment.timeToPass = timeToPass;
segment.oncomingTraffic = oncomingTraffic;
segment.trafficInFront = trafficInFront;

save(strcat(pathToMat(1:end-4), "_withTraffic.mat"), '-struct', 'segment');







