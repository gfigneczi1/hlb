function evaluator_roadGeneration(segments,config)
% This function takes a real-vehicle measurement which contains accurate
% GPS localization data and converts it to a CarMaker road file. 
% The basis of the conversion is the mid-lane UTF coordinates.

% Segmentor: segmentor_driverModel

% @copyright (C) 2022 Robert Bosch GmbH.
% The reproduction, distribution and utilization of this file as
% well as the communication of its contents to others without express
% authorization is prohibited. Offenders will be held liable for the
% payment of damages. All rights reserved in the event of the grant
% of a patent, utility model or design.
% @version 1.0
% @date 2022-12-07

cuttingStart = 11300;
cuttingEnd = 12150;

segment = segments.uncut;

segment.X_abs = segment.LongPos_abs * 40075000 .* cos(segment.LatPos_abs*pi()/180) / 360;
segment.Y_abs = segment.LatPos_abs * 111.32*1000;
segment.theta_calc = theta_recalc(segment)';

%% Get the mid lane data during the whole measurement
for j=1:size(segment.c2,1)    
    corridor(j,1:2) = pos_tf2GPS(segment.X_abs(j),segment.Y_abs(j),segment.theta_calc(j),0.5*(segment.c01_left(j)+segment.c01_right(j)));
end

%% Generate Road File
% Cutting
cut(:,1:2) = corridor(cuttingStart:cuttingEnd,1:2);
corridor = cut;

% Preparing RoadNetworkLength
length = 0;
for j=2:size(corridor,1)
    length = length + ((corridor(j,1)-corridor(j-1,1))^2+(corridor(j,2)-corridor(j-1,2))^2)^(1/2);
end   

% Preparing BBox
minZ = -11;
maxZ = 11;
minX = corridor(1,1);
maxX = corridor(2,1);
minY = corridor(1,2);
maxY = corridor(2,2);
for j=2:size(corridor,1)
    if corridor(j,1) < minX
        minX = corridor(j,1);
    end
    if corridor(j,1) > maxX
        maxX = corridor(j,1);
    end
    if corridor(j,2) < minY
        minY = corridor(j,2);
    end
    if corridor(j,2) > maxY
        maxY = corridor(j,2);
    end
end

% RoadFile
my_dir  = pwd;
road_dir = my_dir + "\04_simulation\vmcsim.platform\Data\Road\roadfile.rd5";
road = fopen(road_dir, 'w' );
    fprintf(road, '#INFOFILE1.1 - Do not remove this line!\n');
    fprintf(road, 'FileIdent = IPGRoad 8.0\n');
    fprintf(road, 'FileCreator = CarMaker 8.1.1\n');
    fprintf(road, 'LibVersion = 8.1.1\n');
    fprintf(road, 'Country = DEU\n');
    fprintf(road, 'nLinks = 1\n');
    fprintf(road, 'nJunctions = 0\n');
    fprintf(road, 'nObjects = \n');
    fprintf(road, 'nRoutes = 0\n');
    fprintf(road, 'RoadNetworkLength = %f\n', length);
    fprintf(road, 'BBox = %d %d %d %d %d %d\n', minX, maxX, minY, maxY, minZ, maxZ);
    fprintf(road, 'RST.Unit = kmh\n');
    fprintf(road, 'RST = -1 -1 -1 -1 -1 -1 -1 -1\n'); %speed limit (negative -> unlimited)
    fprintf(road, 'Movie = 0.2 1 0.01 1.5 1.5 1 1\n'); %visualisation parameters
    fprintf(road, 'PathMode = -1\n'); %Driving side of the road (-1 -> right; 1 -> left)
    fprintf(road, 'Link.0.ID = 0\n');
    fprintf(road, 'Link.0.Junctions = -1 -1 -2 -1\n');
    phi = 0;
    if ((corridor(2,2)-corridor(1,2))>0 && ((corridor(2,1)-corridor(1,1))<0))
        phi = 180;  %for the arctan of the Start Angle
    end
    fprintf(road, 'Link.0.Node0 = %f %f 0 %f\n', 0, 0, atand((corridor(2,2)-corridor(1,2))/(corridor(2,1)-corridor(1,1)))+phi);
    fprintf(road, 'Link.0.RL.ID = 1\n');
    fprintf(road, 'Link.0.Seg.0.ID = 5\n');
    fprintf(road, 'Link.0.Seg.0.Type = PointList\n');
    
    smoothing = 1; %0.....1 (0 <-> no smoothing)
    fprintf(road, 'Link.0.Seg.0.Param = %f %f 0 %f 0 0 0 0\n' , corridor(size(corridor,1)-1,2) - corridor(size(corridor,1),2), corridor(size(corridor,1)-1,1) - corridor(size(corridor,1),1), smoothing);
    fprintf(road, 'Link.0.Seg.0.PointList:\n');
    
    %Generate Road Points
    step = 20;
    for j=1:step:size(corridor,1)
        fprintf(road, '\t%f %f\n',corridor(j,1)-corridor(1,1),corridor(j,2)-corridor(1,2));
    end  
    
    %RoadType
    widthLane = 3.5;
    ID = 10;
    start = 0;
    if segment.roadType(2) == 1
        fprintf(road, 'Link.0.LaneSection.0.ID = %i\n', ID);
        ID = ID + 1;
        fprintf(road, 'Link.0.LaneSection.0.Start = %f\n', start);
        fprintf(road, 'Link.0.LaneSection.0.LaneL.0.ID = %i\n', ID);
        LaneL0_ID = ID;
        ID = ID + 1;
        fprintf(road, 'Link.0.LaneSection.0.LaneL.0 = 0 %f %f 0 0 0 0\n', widthLane, widthLane);
        fprintf(road, 'Link.0.LaneSection.0.LaneL.0.ARP = %i %i %i %i %i %i\n', ID, ID+1, ID+2, ID+3, ID+4, ID+5);
        ID = ID + 6;
        fprintf(road, 'Link.0.LaneSection.0.LaneL.1.ID = %i\n', ID);
        LaneL1_ID = ID;
        ID = ID + 1;
        fprintf(road, 'Link.0.LaneSection.0.LaneL.1 = 0 %f %f 4 0 0 0\n', widthLane, widthLane);
        fprintf(road, 'Link.0.LaneSection.0.LaneL.1.ARP = %i %i %i %i %i %i\n', ID, ID+1, ID+2, ID+3, ID+4, ID+5);
        ID = ID + 6;
        fprintf(road, 'Link.0.LaneSection.0.LaneL.2.ID = %i\n', ID);
        ID = ID + 1;
        fprintf(road, 'Link.0.LaneSection.0.LaneL.2 = 0 2.5 2.5 5 0 0 0\n');
        fprintf(road, 'Link.0.LaneSection.0.LaneR.0.ID = %i\n', ID);
        LaneR0_ID = ID;
        ID = ID + 1;
        fprintf(road, 'Link.0.LaneSection.0.LaneR.0 = 0 %f %f 0 0 0 0\n', widthLane, widthLane);
        fprintf(road, 'Link.0.LaneSection.0.LaneR.0.ARP = %i %i %i %i %i %i\n', ID, ID+1, ID+2, ID+3, ID+4, ID+5);
        ID = ID + 6;
        fprintf(road, 'Link.0.LaneSection.0.LaneR.1.ID = %i\n', ID);
        LaneR1_ID = ID;
        ID = ID + 1;
        fprintf(road, 'Link.0.LaneSection.0.LaneR.1 = 0 %f %f 4 0 0 0\n', widthLane, widthLane);
        fprintf(road, 'Link.0.LaneSection.0.LaneR.1.ARP = %i %i %i %i %i %i\n', ID, ID+1, ID+2, ID+3, ID+4, ID+5);
        ID = ID + 6;
        fprintf(road, 'Link.0.LaneSection.0.LaneR.2.ID = %i\n', ID);
        ID = ID + 1;
        fprintf(road, 'Link.0.LaneSection.0.LaneR.2 = 0 2.5 2.5 5 0 0 0\n');
        
        fprintf(road, 'LanePath.0 = %i %i 0.25 10 0.1 0.1\n', ID, LaneL0_ID);
        fprintf(road, 'LanePath.1 = %i %i 2 10 0.1 0.1\n', ID+1, LaneL1_ID);
        fprintf(road, 'LanePath.2 = %i %i 2 10 0.1 0.1\n', ID+2, LaneR0_ID);
        fprintf(road, 'LanePath.3 = %i %i 0.25 10 0.1 0.1\n', ID+3, LaneR1_ID);
        ID = ID + 4;
        RouteID = ID - 2;
        
        fprintf(road, 'Route.0.ID = %i\n', ID);
        ID = ID + 1;
        fprintf(road, 'Route.0.Name = Route_0\n');
        fprintf(road, 'Route.0.DrvPath.ID = %i\n', ID);
        ID = ID + 1;
        fprintf(road, 'Route.0.DrvPath: \n\t%i\n', RouteID);
    
        fprintf(road, 'RL.1.RoadMarking.0.ID = %i %i\n', ID, LaneL0_ID); % gestichelt (Mitte)
        ID = ID + 1;
        fprintf(road, 'RL.1.RoadMarking.0 = 0 0 0 1 0 -1 0.15 0 2 0 0 2 4 1 1 0 ""\n');
        fprintf(road, 'RL.1.RoadMarking.0.Material.0 = 1.0,1.0,1.0 0 0 0 0 0 0 0 0 0 0 0\n');

        fprintf(road, 'RL.1.RoadMarking.1.ID = %i %i\n', ID, LaneL1_ID); % durchgezogen
        ID = ID + 1;
        fprintf(road, 'RL.1.RoadMarking.1 = 0 0 0 1 0 -1 0.15 0 1 0 0 2 4 1 1 0 ""\n');
        fprintf(road, 'RL.1.RoadMarking.1.Material.0 = 1.0,1.0,1.0 0 0 0 0 0 0 0 0 0 0 0\n');
        
        fprintf(road, 'RL.1.RoadMarking.2.ID = %i %i\n', ID, LaneR0_ID); % durchgezogen
        ID = ID + 1;
        fprintf(road, 'RL.1.RoadMarking.2 = 0 0 0 1 0 -1 0.15 0 1 0 0 2 4 1 1 0 ""\n');
        fprintf(road, 'RL.1.RoadMarking.2.Material.0 = 1.0,1.0,1.0 0 0 0 0 0 0 0 0 0 0 0\n');
        
    elseif segment.roadType(2) == 2
        roadType = 1;
        a = 0;
        b = 0;
        for j=1:size(segment.c02_right)
            if (segment.c02_left(j) ~= 0 && a == 0)
                roadType = roadType + 1;
                a = 1;
            elseif (segment.c02_left(j) ~= 0 && b == 0)
                roadType = roadType + 1;
                b = 1;
            end   
        end
        if roadType == 2
            fprintf(road, 'Link.0.LaneSection.0.ID = %i\n', ID);
            ID = ID + 1;
            fprintf(road, 'Link.0.LaneSection.0.Start = %f\n', start);
            fprintf(road, 'Link.0.LaneSection.0.LaneL.0.ID = %i\n', ID);
            LaneL0_ID = ID;
            ID = ID + 1;
            fprintf(road, 'Link.0.LaneSection.0.LaneL.0 = 0 %f %f 0 0 0 0\n', widthLane, widthLane);
            fprintf(road, 'Link.0.LaneSection.0.LaneL.0.ARP = %i %i %i %i %i %i\n', ID, ID+1, ID+2, ID+3, ID+4, ID+5);
            ID = ID + 6;
            fprintf(road, 'Link.0.LaneSection.0.LaneL.1.ID = %i\n', ID);
            LaneL1_ID = ID;
            ID = ID + 1;
            fprintf(road, 'Link.0.LaneSection.0.LaneL.1 = 0 %f %f 0 0 0 0\n', widthLane, widthLane);
            fprintf(road, 'Link.0.LaneSection.0.LaneL.1.ARP = %i %i %i %i %i %i\n', ID, ID+1, ID+2, ID+3, ID+4, ID+5);
            ID = ID + 6;
            fprintf(road, 'Link.0.LaneSection.0.LaneL.2.ID = %i\n', ID);
            LaneL2_ID = ID;
            ID = ID + 1;
            fprintf(road, 'Link.0.LaneSection.0.LaneL.2 = 0 %f %f 4 0 0 0\n', widthLane, widthLane);
            fprintf(road, 'Link.0.LaneSection.0.LaneL.2.ARP = %i %i %i %i %i %i\n', ID, ID+1, ID+2, ID+3, ID+4, ID+5);
            ID = ID + 6;
            fprintf(road, 'Link.0.LaneSection.0.LaneL.3.ID = %i\n', ID);
            ID = ID + 1;
            fprintf(road, 'Link.0.LaneSection.0.LaneL.3 = 0 2.5 2.5 5 0 0 0\n');
            fprintf(road, 'Link.0.LaneSection.0.LaneR.0.ID = %i\n', ID);
            LaneR0_ID = ID;
            ID = ID + 1;
            fprintf(road, 'Link.0.LaneSection.0.LaneR.0 = 0 %f %f 0 0 0 0\n', widthLane, widthLane);
            fprintf(road, 'Link.0.LaneSection.0.LaneR.0.ARP = %i %i %i %i %i %i\n', ID, ID+1, ID+2, ID+3, ID+4, ID+5);
            ID = ID + 6;
            fprintf(road, 'Link.0.LaneSection.0.LaneR.1.ID = %i\n', ID);
            LaneR1_ID = ID;
            ID = ID + 1;
            fprintf(road, 'Link.0.LaneSection.0.LaneR.1 = 0 %f %f 0 0 0 0\n', widthLane, widthLane);
            fprintf(road, 'Link.0.LaneSection.0.LaneR.1.ARP = %i %i %i %i %i %i\n', ID, ID+1, ID+2, ID+3, ID+4, ID+5);
            ID = ID + 6;
            fprintf(road, 'Link.0.LaneSection.0.LaneR.2.ID = %i\n', ID);
            LaneR2_ID = ID;
            ID = ID + 1;
            fprintf(road, 'Link.0.LaneSection.0.LaneR.2 = 0 %f %f 4 0 0 0\n', widthLane, widthLane);
            fprintf(road, 'Link.0.LaneSection.0.LaneR.2.ARP = %i %i %i %i %i %i\n', ID, ID+1, ID+2, ID+3, ID+4, ID+5);
            ID = ID + 6;
            fprintf(road, 'Link.0.LaneSection.0.LaneR.3.ID = %i\n', ID);
            ID = ID + 1;
            fprintf(road, 'Link.0.LaneSection.0.LaneR.3 = 0 2.5 2.5 5 0 0 0\n');

            fprintf(road, 'LanePath.0 = %i %i 0.25 10 0.1 0.1\n', ID, LaneL0_ID);
            fprintf(road, 'LanePath.1 = %i %i 2 10 0.1 0.1\n', ID+1, LaneL1_ID);
            fprintf(road, 'LanePath.2 = %i %i 2 10 0.1 0.1\n', ID+2, LaneL2_ID);
            fprintf(road, 'LanePath.3 = %i %i 2 10 0.1 0.1\n', ID+3, LaneR0_ID);
            fprintf(road, 'LanePath.4 = %i %i 2 10 0.1 0.1\n', ID+4, LaneR1_ID);
            fprintf(road, 'LanePath.5 = %i %i 0.25 10 0.1 0.1\n', ID+5, LaneR2_ID);
            ID = ID + 6;
            RouteID = ID - 2;
            
            fprintf(road, 'Route.0.ID = %i\n', ID);
            ID = ID + 1;
            fprintf(road, 'Route.0.Name = Route_0\n');
            fprintf(road, 'Route.0.DrvPath.ID = %i\n', ID);
            ID = ID + 1;
            fprintf(road, 'Route.0.DrvPath: \n\t%i\n', RouteID);
    
            fprintf(road, 'RL.1.RoadMarking.0.ID = %i %i\n', ID, LaneL0_ID); % doppelt (Mitte)
            ID = ID + 1;
            fprintf(road, 'RL.1.RoadMarking.0 = 0 0 0 1 0 -1 0.15 0 4 0 0 2 4 1 1 0 ""\n');
            fprintf(road, 'RL.1.RoadMarking.0.Material.0 = 1.0,1.0,1.0 0 0 0 0 0 0 0 0 0 0 0\n');
            fprintf(road, 'RL.1.RoadMarking.1.ID = %i %i\n', ID, LaneL1_ID); % gestichelt 
            ID = ID + 1;
            fprintf(road, 'RL.1.RoadMarking.1 = 0 0 0 1 0 -1 0.15 0 2 0 0 2 4 1 1 0 ""\n');
            fprintf(road, 'RL.1.RoadMarking.1.Material.0 = 1.0,1.0,1.0 0 0 0 0 0 0 0 0 0 0 0\n');
            fprintf(road, 'RL.1.RoadMarking.2.ID = %i %i\n', ID, LaneL2_ID); % durchgezogen
            ID = ID + 1;
            fprintf(road, 'RL.1.RoadMarking.2 = 0 0 0 1 0 -1 0.15 0 1 0 0 2 4 1 1 0 ""\n');
            fprintf(road, 'RL.1.RoadMarking.2.Material.0 = 1.0,1.0,1.0 0 0 0 0 0 0 0 0 0 0 0\n');
            fprintf(road, 'RL.1.RoadMarking.3.ID = %i %i\n', ID, LaneR0_ID); % gestichelt 
            ID = ID + 1;
            fprintf(road, 'RL.1.RoadMarking.3 = 0 0 0 1 0 -1 0.15 0 2 0 0 2 4 1 1 0 ""\n');
            fprintf(road, 'RL.1.RoadMarking.3.Material.0 = 1.0,1.0,1.0 0 0 0 0 0 0 0 0 0 0 0\n');
            fprintf(road, 'RL.1.RoadMarking.4.ID = %i %i\n', ID, LaneR1_ID); % durchgezogen
            ID = ID + 1;
            fprintf(road, 'RL.1.RoadMarking.4 = 0 0 0 1 0 -1 0.15 0 1 0 0 2 4 1 1 0 ""\n');
            fprintf(road, 'RL.1.RoadMarking.4.Material.0 = 1.0,1.0,1.0 0 0 0 0 0 0 0 0 0 0 0\n');

        elseif roadType == 3
            fprintf(road, 'Link.0.LaneSection.0.ID = %i\n', ID);
            ID = ID + 1;
            fprintf(road, 'Link.0.LaneSection.0.Start = %f\n', start);
            fprintf(road, 'Link.0.LaneSection.0.LaneL.0.ID = %i\n', ID);
            LaneL0_ID = ID;
            ID = ID + 1;
            fprintf(road, 'Link.0.LaneSection.0.LaneL.0 = 0 %f %f 0 0 0 0\n', widthLane, widthLane);
            fprintf(road, 'Link.0.LaneSection.0.LaneL.0.ARP = %i %i %i %i %i %i\n', ID, ID+1, ID+2, ID+3, ID+4, ID+5);
            ID = ID + 6;
            fprintf(road, 'Link.0.LaneSection.0.LaneL.1.ID = %i\n', ID);
            LaneL1_ID = ID;
            ID = ID + 1;
            fprintf(road, 'Link.0.LaneSection.0.LaneL.1 = 0 %f %f 0 0 0 0\n', widthLane, widthLane);
            fprintf(road, 'Link.0.LaneSection.0.LaneL.1.ARP = %i %i %i %i %i %i\n', ID, ID+1, ID+2, ID+3, ID+4, ID+5);
            ID = ID + 6;
            fprintf(road, 'Link.0.LaneSection.0.LaneL.2.ID = %i\n', ID);
            LaneL2_ID = ID;
            ID = ID + 1;
            fprintf(road, 'Link.0.LaneSection.0.LaneL.2 = 0 %f %f 0 0 0 0\n', widthLane, widthLane);
            fprintf(road, 'Link.0.LaneSection.0.LaneL.2.ARP = %i %i %i %i %i %i\n', ID, ID+1, ID+2, ID+3, ID+4, ID+5);
            ID = ID + 6;
            fprintf(road, 'Link.0.LaneSection.0.LaneL.3.ID = %i\n', ID);
            LaneL3_ID = ID;
            ID = ID + 1;
            fprintf(road, 'Link.0.LaneSection.0.LaneL.3 = 0 %f %f 4 0 0 0\n', widthLane, widthLane);
            fprintf(road, 'Link.0.LaneSection.0.LaneL.3.ARP = %i %i %i %i %i %i\n', ID, ID+1, ID+2, ID+3, ID+4, ID+5);
            ID = ID + 6;
            fprintf(road, 'Link.0.LaneSection.0.LaneL.4.ID = %i\n', ID);
            ID = ID + 1;
            fprintf(road, 'Link.0.LaneSection.0.LaneL.4 = 0 2.5 2.5 5 0 0 0\n');
            fprintf(road, 'Link.0.LaneSection.0.LaneR.0.ID = %i\n', ID);
            LaneR0_ID = ID;
            ID = ID + 1;
            fprintf(road, 'Link.0.LaneSection.0.LaneR.0 = 0 %f %f 0 0 0 0\n', widthLane, widthLane);
            fprintf(road, 'Link.0.LaneSection.0.LaneR.0.ARP = %i %i %i %i %i %i\n', ID, ID+1, ID+2, ID+3, ID+4, ID+5);
            ID = ID + 6;
            fprintf(road, 'Link.0.LaneSection.0.LaneR.1.ID = %i\n', ID);
            LaneR1_ID = ID;
            ID = ID + 1;
            fprintf(road, 'Link.0.LaneSection.0.LaneR.1 = 0 %f %f 0 0 0 0\n', widthLane, widthLane);
            fprintf(road, 'Link.0.LaneSection.0.LaneR.1.ARP = %i %i %i %i %i %i\n', ID, ID+1, ID+2, ID+3, ID+4, ID+5);
            ID = ID + 6;
            fprintf(road, 'Link.0.LaneSection.0.LaneR.2.ID = %i\n', ID);
            LaneR2_ID = ID;
            ID = ID + 1;
            fprintf(road, 'Link.0.LaneSection.0.LaneR.2 = 0 %f %f 0 0 0 0\n', widthLane, widthLane);
            fprintf(road, 'Link.0.LaneSection.0.LaneR.2.ARP = %i %i %i %i %i %i\n', ID, ID+1, ID+2, ID+3, ID+4, ID+5);
            ID = ID + 6;
            fprintf(road, 'Link.0.LaneSection.0.LaneR.3.ID = %i\n', ID);
            LaneR3_ID = ID;
            ID = ID + 1;
            fprintf(road, 'Link.0.LaneSection.0.LaneR.3 = 0 %f %f 4 0 0 0\n', widthLane, widthLane);
            fprintf(road, 'Link.0.LaneSection.0.LaneR.3.ARP = %i %i %i %i %i %i\n', ID, ID+1, ID+2, ID+3, ID+4, ID+5);
            ID = ID + 6;
            fprintf(road, 'Link.0.LaneSection.0.LaneR.4.ID = %i\n', ID);
            ID = ID + 1;
            fprintf(road, 'Link.0.LaneSection.0.LaneR.4 = 0 2.5 2.5 5 0 0 0\n');

            fprintf(road, 'LanePath.0 = %i %i 0.25 10 0.1 0.1\n', ID, LaneL0_ID);
            fprintf(road, 'LanePath.1 = %i %i 2 10 0.1 0.1\n', ID+1, LaneL1_ID);
            fprintf(road, 'LanePath.2 = %i %i 2 10 0.1 0.1\n', ID+2, LaneL2_ID);
            fprintf(road, 'LanePath.3 = %i %i 2 10 0.1 0.1\n', ID+3, LaneL3_ID);
            fprintf(road, 'LanePath.4 = %i %i 2 10 0.1 0.1\n', ID+4, LaneR0_ID);
            fprintf(road, 'LanePath.5 = %i %i 2 10 0.1 0.1\n', ID+5, LaneR1_ID);
            fprintf(road, 'LanePath.6 = %i %i 2 10 0.1 0.1\n', ID+6, LaneR2_ID);
            fprintf(road, 'LanePath.7 = %i %i 0.25 10 0.1 0.1\n', ID+7, LaneR3_ID);
            ID = ID + 8;
            RouteID = ID - 3;
            
            fprintf(road, 'Route.0.ID = %i\n', ID);
            ID = ID + 1;
            fprintf(road, 'Route.0.Name = Route_0\n');
            fprintf(road, 'Route.0.DrvPath.ID = %i\n', ID);
            ID = ID + 1;
            fprintf(road, 'Route.0.DrvPath: \n\t%i\n', RouteID);
    
            fprintf(road, 'RL.1.RoadMarking.0.ID = %i %i\n', ID, LaneL0_ID); % doppelt (Mitte)
            ID = ID + 1;
            fprintf(road, 'RL.1.RoadMarking.0 = 0 0 0 1 0 -1 0.15 0 4 0 0 2 4 1 1 0 ""\n');
            fprintf(road, 'RL.1.RoadMarking.0.Material.0 = 1.0,1.0,1.0 0 0 0 0 0 0 0 0 0 0 0\n');
            fprintf(road, 'RL.1.RoadMarking.1.ID = %i %i\n', ID, LaneL1_ID); % gestichelt 
            ID = ID + 1;
            fprintf(road, 'RL.1.RoadMarking.1 = 0 0 0 1 0 -1 0.15 0 2 0 0 2 4 1 1 0 ""\n');
            fprintf(road, 'RL.1.RoadMarking.1.Material.0 = 1.0,1.0,1.0 0 0 0 0 0 0 0 0 0 0 0\n');
            fprintf(road, 'RL.1.RoadMarking.2.ID = %i %i\n', ID, LaneL2_ID); % gestichelt
            ID = ID + 1;
            fprintf(road, 'RL.1.RoadMarking.2 = 0 0 0 1 0 -1 0.15 0 2 0 0 2 4 1 1 0 ""\n');
            fprintf(road, 'RL.1.RoadMarking.2.Material.0 = 1.0,1.0,1.0 0 0 0 0 0 0 0 0 0 0 0\n');
            fprintf(road, 'RL.1.RoadMarking.3.ID = %i %i\n', ID, LaneL3_ID); % durchgezogen
            ID = ID + 1;
            fprintf(road, 'RL.1.RoadMarking.3 = 0 0 0 1 0 -1 0.15 0 1 0 0 2 4 1 1 0 ""\n');
            fprintf(road, 'RL.1.RoadMarking.3.Material.0 = 1.0,1.0,1.0 0 0 0 0 0 0 0 0 0 0 0\n');
            fprintf(road, 'RL.1.RoadMarking.4.ID = %i %i\n', ID, LaneR0_ID); % gestichelt 
            ID = ID + 1;
            fprintf(road, 'RL.1.RoadMarking.4 = 0 0 0 1 0 -1 0.15 0 2 0 0 2 4 1 1 0 ""\n');
            fprintf(road, 'RL.1.RoadMarking.4.Material.0 = 1.0,1.0,1.0 0 0 0 0 0 0 0 0 0 0 0\n');
            fprintf(road, 'RL.1.RoadMarking.5.ID = %i %i\n', ID, LaneR1_ID); % gestichelt
            ID = ID + 1;
            fprintf(road, 'RL.1.RoadMarking.5 = 0 0 0 1 0 -1 0.15 0 2 0 0 2 4 1 1 0 ""\n');
            fprintf(road, 'RL.1.RoadMarking.5.Material.0 = 1.0,1.0,1.0 0 0 0 0 0 0 0 0 0 0 0\n');
            fprintf(road, 'RL.1.RoadMarking.6.ID = %i %i\n', ID, LaneR2_ID); % durchgezogen
            ID = ID + 1;
            fprintf(road, 'RL.1.RoadMarking.6 = 0 0 0 1 0 -1 0.15 0 1 0 0 2 4 1 1 0 ""\n');
            fprintf(road, 'RL.1.RoadMarking.6.Material.0 = 1.0,1.0,1.0 0 0 0 0 0 0 0 0 0 0 0\n');
        end
    end
    
    fprintf(road, 'MaxUsedObjId = %i\n', ID);
    fclose(road);


end
