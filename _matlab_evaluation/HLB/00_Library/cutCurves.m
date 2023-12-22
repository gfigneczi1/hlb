function [curves] = cutCurves(driverModelOutput)
% based on GPS coordinates
% defining rectangles of curves
% cell array. 1st column: start rectangle, 2nd column: finish rectangle
r{1,1} = [1431760 1431770 5233075 5233085];
r{1,2} = [1430390 1430410 5234144 5234154];
r{2,1} = [1429175 1429185 5235655 5235665];
r{2,2} = [1428698 1428708 5236085 5236095];
r{3,1} = [1428715 1428725 5236075 5236085];
r{3,2} = [1428015 1428025 5236600 5236610];
r{4,1} = [1426155 1426165 5238595 5238605];
r{4,2} = [1423325 1423335 5239395 5239405];
r{5,1} = [1423325 1423335 5239395 5239405];
r{5,2} = [1420995 1421005 5239710 5239720];
r{6,1} = [1420168 1420178 5240336 5240346];
r{6,2} = [1418885 1418895 5240792 5240802];
r{7,1} = [1416945 1416955 5242025 5242035];
r{7,2} = [1414200 4141210 5242968 5242978];
r{8,1} = [1412995 1413005 5243425 5243435];
r{8,2} = [1410615 1410625 5243390 5243400];
r{9,1} = [1409595 1409605 5242815 5242825];
r{9,2} = [1407235 1407245 5243245 5243255];
r{10,1} = [1405665 1405675 5244505 5244495];
r{10,2} = [1403995 1404005 5247145 5247155];

for i=1:size(r,1)
    for j=1:2
        temp = r{i,j};
        r{i,j}(3) = temp(2);
        r{i,j}(4) = temp(1);
        r{i,j}(5:6) = temp(3);
        r{i,j}(7:8) = temp(4);
    end
end

r = extendCurves(r, 50);
vxGlobal = mean(diff(driverModelOutput.traj(:,1)));


for i=1:size(r,1)
    if (vxGlobal > 0)
        starts = inpolygon(driverModelOutput.traj(:,1), driverModelOutput.traj(:,2), ...
            r{i,2}(1:4), r{i,2}(5:8));
        stops = inpolygon(driverModelOutput.traj(:,1), driverModelOutput.traj(:,2), ...
            r{i,1}(1:4), r{i,1}(5:8));
    else
        starts = inpolygon(driverModelOutput.traj(:,1), driverModelOutput.traj(:,2), ...
            r{i,1}(1:4), r{i,1}(5:8));
        stops = inpolygon(driverModelOutput.traj(:,1), driverModelOutput.traj(:,2), ...
            r{i,2}(1:4), r{i,2}(5:8));
    end
    start = find(starts==1,1);
    stop = find(stops==1,1);
    if (isempty(start) || isempty(stop))
        start = 1; stop = size(driverModelOutput.traj,1);
        disp(strcat('no curve found for r(',num2str(i),')'));
    end
    if (vxGlobal > 0)
        curves(size(r,1)-i+1).traj = driverModelOutput.traj(start:stop,:);
        curves(size(r,1)-i+1).ref = driverModelOutput.ref(start:stop,:);
        curves(size(r,1)-i+1).cor = driverModelOutput.cor(start:stop,:);
        curves(size(r,1)-i+1).orient = driverModelOutput.orient(start:stop,:);
        curves(size(r,1)-i+1).corLeft = [driverModelOutput.cor(start:stop,1)-sin(curves(size(r,1)-i+1).orient)*1.875 ...
            driverModelOutput.cor(start:stop,2)+cos(curves(size(r,1)-i+1).orient)*1.875];
        curves(size(r,1)-i+1).corRight = [driverModelOutput.cor(start:stop,1)+sin(curves(size(r,1)-i+1).orient)*1.875 ...
            driverModelOutput.cor(start:stop,2)-cos(curves(size(r,1)-i+1).orient)*1.875];
        curves(size(r,1)-i+1).curv = driverModelOutput.curv(start:stop,:);
        curves(size(r,1)-i+1).displacement = [0; cumtrapz((diff(curves(size(r,1)-i+1).cor(:,1)).^2+diff(curves(size(r,1)-i+1).cor(:,2)).^2).^0.5)];
        %curves(i).replan = driverModelOutput.replan(start:stop);
        if(isfield(driverModelOutput,'GT_U'))
            curves(size(r,1)-i+1).GT_U = driverModelOutput.GT_U(:,start:stop);
        end
        %curves(i).GT_U = curves(i).GT_U(:, curves(i).replan==1);
        if(isfield(driverModelOutput,'GT_Y'))
            curves(size(r,1)-i+1).GT_Y = driverModelOutput.GT_Y(:,start:stop);
        end
        %curves(i).GT_Y = curves(i).GT_Y(:, curves(i).replan==1);
        curves(size(r,1)-i+1).start = start;
        curves(size(r,1)-i+1).stop = stop;
    else
        curves(i).traj = driverModelOutput.traj(start:stop,:);
        curves(i).ref = driverModelOutput.ref(start:stop,:);
        curves(i).cor = driverModelOutput.cor(start:stop,:);
        curves(i).orient = driverModelOutput.orient(start:stop,:);
        curves(i).corLeft = [driverModelOutput.cor(start:stop,1)-sin(curves(i).orient)*1.875 ...
            driverModelOutput.cor(start:stop,2)+cos(curves(i).orient)*1.875];
        curves(i).corRight = [driverModelOutput.cor(start:stop,1)+sin(curves(i).orient)*1.875 ...
            driverModelOutput.cor(start:stop,2)-cos(curves(i).orient)*1.875];
        curves(i).curv = driverModelOutput.curv(start:stop,:);
        curves(i).displacement = [0; cumtrapz((diff(curves(i).cor(:,1)).^2+diff(curves(i).cor(:,2)).^2).^0.5)];
        %curves(i).replan = driverModelOutput.replan(start:stop);
        if(isfield(driverModelOutput,'GT_U'))
            curves(i).GT_U = driverModelOutput.GT_U(:,start:stop);
        end
        %curves(i).GT_U = curves(i).GT_U(:, curves(i).replan==1);
        if(isfield(driverModelOutput,'GT_Y'))
            curves(i).GT_Y = driverModelOutput.GT_Y(:,start:stop);
        end
        %curves(i).GT_Y = curves(i).GT_Y(:, curves(i).replan==1);
        curves(i).start = start;
        curves(i).stop = stop;
    end
end
end

function r = extendCurves(r, extension)
    for i=1:size(r,1)
        % r is CCW
        for j=1:2
            r{i,j}(1) = r{i,j}(1) - extension;
            r{i,j}(4) = r{i,j}(4) - extension;
            r{i,j}(2:3) = r{i,j}(2:3) + extension;
            r{i,j}(5:6) = r{i,j}(5:6) - extension;
            r{i,j}(7:8) = r{i,j}(7:8) + extension;
        end
    end
end
