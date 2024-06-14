function [traj, ref, cor, corLeft, corRight, orient, curv, GT_U, GT_Y, validPoints, replan_array, GT_U_array, GT_Y_array, trajFull, refFull, corFull] = functional_driverModelWithPlanner(x, use)
    % This function aims to calculate the trajectory over a route (e.g.,
    % many kilometres) based on the Linear Driver Model (LDM).
    % This model is discovered by Gergo Ferenc Igneczi. It is a
    % parametrizable model through the coefficients of the linear
    % functions. The model uses the corridor geometry in front of the
    % vehicle and calculates the node point offsets. As it returns only
    % node points, a trajectory planner must be used which fits a curve on 
    % the node points. In this configuration 3 node points are used. The
    % fourth (the zero one) is always in the planner frame origo. 
    % Planner frame: the frame the model requests the corridor in and
    % returns the node points in. It can be any 2D frame described by a 3
    % element pose (X,Y, theta).
    % All orientation is understood as angles in radian.
    % Current selection for the planner frame: it is positioned within the
    % Global UTF frame, always selected to produce continuous global
    % trajectories (upon replanning, it is positioned to the last valid
    % trajectory point).
    % Usage:
    % - from a global variable the startIdx and endIdx must be given. These
    % are the indices which denote the section of a complete measurement to
    % be evaluated. These can cover hundreds of kilometres.
    % - from a global variable the lookAheadDistance (in meters) must be
    % given. Based on the lookAheadDistance and the actual point in the
    % measurement endIndex is calculated, which is the end of the
    % subsegment within the measurement section. If endIndex is bigger than
    % the number of samples in the measurement, the iterations are stopped,
    % while if endIndex > endIdx the last trajectory will be completely
    % saved and stored.
    
    % Inputs:
    % - P: 21*1 parameter vector containing the linear coefficients of the
    % LDM
    % - use: how the planner is used together with the model: use = 0:
    % ground-truth from the measurement is taken for the node point
    % calculation (= curve fitting on the measured data), use = 1: LDM is
    % used to calculate the node points.
    % Outputs:
    % - traj: planned trajectory in the global frame
    % - ref: measured reference trajectory = trajectory driven by human
    % drivers
    % - cor: global path of the midlane
    % - corLeft/corRight: left and right borders degined globally
    % - orient: orientation of the trajectory in the ego frame
    % - curv: curvature of the planned trajectory in the ego frame
    % - GT_U: array of the input variables, 7xn, 7: 3 average curvature and
    % 4 curvature gradients in the node points, n = number of trajectory
    % replanning
    % - GT_Y: similarly to GT_U, array of the node point characteristics
    
    % @copyright (C) 2022 Robert Bosch GmbH.
    % The reproduction, distribution and utilization of this file as
    % well as the communication of its contents to others without express
    % authorization is prohibited. Offenders will be held liable for the
    % payment of damages. All rights reserved in the event of the grant
    % of a patent, utility model or design.
    % @version 1.0
    % @date 2022-12-07
    
    clearvars GT_U GT_Y GT_theta anchorPoints_array 
    global segment_m startIdx lookAheadDistance endIdx indexes P3in model_ver
    
    %% Initialization
    P = x;
    traj = []; % empty array
    ref = [];
    cor = [];
    curv = [];
    orient = [];
    kappa_array = [];
    GT_Y = [];
    GT_theta = [];
    GT_U = [];
    anchorPoints_array = [];
    n = size(segment_m,1);
    endIdx = min(n,endIdx);
    startingIndex = startIdx;
    endIndex =  1000;
    replanCycle = 10;
    replan = 1;
    initialized = 0;
    counter = 1; % counter, how far we went in the previous trajectory
    validPoints = zeros(size(segment_m,1),1);
    trajFull = [];
    corFull = [];
    refFull = [];
    replan_array = [];
    GT_U_array = [];
    GT_Y_array = [];
    
    %% Iteration through measurement data
    maxStepSize = segment_m(:,indexes.VelocityX) * 0.1; % v * t, t = 100ms (50ms sampling time is assumed);
    for i=startIdx:endIdx
        %replan_array(i-startIdx+1) = replan;
        if ((((segment_m(i,indexes.X_abs)-segment_m(i-1,indexes.X_abs))^2+(segment_m(i,indexes.Y_abs)-segment_m(i-1,indexes.Y_abs))^2)^0.5) > maxStepSize(i))
            % regardless where we are in the planning sequence we trigger a
            % reinit, as a jump in the measurement has been detected
            replan = 1;
            initialized = 0;
            counter = 1; % counter, how far we went in the previous trajectory
        end
        if (initialized && size(trajectoryPreviousPlannerFrame,1) < 0.5*lookAheadDistance/maxStepSize(i))
            replan = 1;
            initialized = 0;
            counter = 1; % counter, how far we went in the previous trajectory
        end

        % if branch: first iteration or replan trigger
        if (replan > 0 || ~initialized)
            if (initialized)
                plannerFrame = getPlannerFrame(segmentsParam, previousPlannerFrame, trajectoryPreviousPlannerFrame, counter, anchorPoints, indeces);
                %plannerFrame = correctPlannerFrame(segment_m(:,indexes.X_abs), segment_m(:,indexes.Y_abs), plannerFrame);
            else
                plannerFrame = [segment_m(i,indexes.X_abs); segment_m(i,indexes.Y_abs); segment_m(i,indexes.theta_calc)];
            end
            previousPlannerFrame = plannerFrame; % returning the last planner frame
            initialized = 1;

            distances = ((segment_m(1:end, indexes.corrX)-segment_m(i,indexes.X_abs)).^2+(segment_m(1:end, indexes.corrY)-segment_m(i,indexes.Y_abs)).^2).^0.5;
            startingIndex = find(distances == min(distances),1);
            distances = distances(startingIndex:end);
            endIndex = find(distances > lookAheadDistance,1);
            
            if (~isempty(endIndex))
                if (startingIndex + endIndex(1,1)) > size(segment_m,1)
                    % hard check: not enough sample points left
                    break;
                end
            end
            
            if (~isempty(endIndex) && ~isempty(startingIndex) && min(distances) <= 40)
                if (distances(endIndex) < 1.1*lookAheadDistance && ...
                        distances(endIndex) > 0.5 * lookAheadDistance && ...
                        ~any(abs(diff(segment_m(startingIndex:(startingIndex+endIndex(1,1)),indexes.X_abs)))>10) && ...
                        ~any(abs(diff(segment_m(startingIndex:(startingIndex+endIndex(1,1)),indexes.Y_abs)))>10))
                    % 0. resetting counter
                    counter = 1;
                    % 1. cutting the subsegment
                    % subsegment = cutSubsegment(segment, startingIndex, startingIndex+endIndex(1,1));
                    subsegment_m = segment_m(startingIndex:(startingIndex+endIndex(1,1)),1:end);
                    % 2. processing the corridor in the k(th) planning
                    % frame
                   
                    % this returns the corridor in the planner frame with
                    % the original sample locations
                    corridorPlannerFrame = corridorPreproc(subsegment_m, plannerFrame);
                    
                    % 3. reference vector calculation
                    % this returns the ego path in the planner frame
                    [refPlannerFrame, refOrientationPlannerFrame] = referenceCalculation(subsegment_m(:,indexes.X_abs), subsegment_m(:,indexes.Y_abs), subsegment_m(:,indexes.theta_calc), plannerFrame);
                    
                    % calculating valid indeces
                    % first points are taken out if those are negative
                    positive_indeces = refPlannerFrame(:,1) >=0;
                    
                    % 4. resample the corridor according to reference X
                    % coordinates
                    corridorPlannerFrameResampled = corridorPlannerFrame; %corridor_resample(corridorPlannerFrame,refPlannerFrame(:,1),subsegment_m(:,18));

                    corridorPlannerFrameResampled = corridorPlannerFrameResampled(positive_indeces,:);
                    refPlannerFrame = refPlannerFrame(positive_indeces,:);
                    refOrientationPlannerFrame = refOrientationPlannerFrame(positive_indeces,:);
                    
                    % 5. calculating the driver model outputs
                    [X, Y, theta, indeces, ~, U] = driverModelPlanner(subsegment_m, P, corridorPlannerFrameResampled, P3in, model_ver);
                    Np = length(indeces);
                    if (any(diff(indeces)==0))
                        % some break in the measurement is found, re-init
                        % everything
                        % TODO: check if this if/else is needed, maybe
                        % duplicated
                        replan = 7;
                        initialized = 0;
                        counter = 1; % counter, how far we went in the previous trajectory
                    else
                        
                        for np=2:length(indeces)
                            % 6. calculation of new trajectory (in local frame)
                            index = find(refPlannerFrame(:,1) >=X(np),1);
                            refPlannerFrame(indeces(np),2) = interp1(refPlannerFrame(index-1:index,1),refPlannerFrame(index-1:index,2),X(np));
                            refPlannerFrame(indeces(np),1) = X(np);
                        end
                        
                        switch use
                            case 0
                                % ground truth
                                Y(2:end) = refPlannerFrame(indeces(2:end),2);
                                theta(2:end) = atan(tan(refOrientationPlannerFrame(indeces(2:end))));
                                anchorPoints = [X; Y; theta];
                                GT_Y = [GT_Y Y-corridorPlannerFrameResampled(indeces,2)];
                                GT_theta = [GT_theta theta];
                                anchorPoints_array = [anchorPoints_array anchorPoints(Np+1:Np+Np)-corridorPlannerFrameResampled(indeces,2)];
                                GT_U = [GT_U U];
                            case 1
                                % LDM
                                anchorPoints = [X; Y; theta];
                                GT_Y = [GT_Y Y-corridorPlannerFrameResampled(indeces,2)];
                                GT_theta = [GT_theta theta];
                                anchorPoints_array = [anchorPoints_array anchorPoints(5:8)-corridorPlannerFrameResampled(indeces,2)];
                                GT_U = [GT_U U];
                            case 2
                                % undisturbed by opposite traffic
                                if (any(segment_m(i:i+indeces(4),indexes.vehicleOppositeSideDisturbing) > 0))
                                    U = zeros(7,1);
                                    Y(2:end) = 0;
                                    theta(2:end) = 0;
                                else
                                    Y(2:end) = refPlannerFrame(indeces(2:end),2);
                                    theta(2:end) = atan(tan(refOrientationPlannerFrame(indeces(2:end))));
                                end
                                anchorPoints = [X; Y; theta];
                                GT_Y = [GT_Y Y-corridorPlannerFrameResampled(indeces,2)];
                                GT_theta = [GT_theta theta];
                                anchorPoints_array = [anchorPoints_array anchorPoints(5:8)-corridorPlannerFrameResampled(indeces,2)];
                                GT_U = [GT_U U];
                            case 3
                                % disturbed by opposite traffic
                                if (any(segment_m(i:i+indeces(4),indexes.vehicleOppositeSideDisturbing) > 0))
                                    Y(2:end) = refPlannerFrame(indeces(2:end),2);
                                    theta(2:end) = atan(tan(refOrientationPlannerFrame(indeces(2:end))));
                                else
                                    U = zeros(7,1);
                                    Y(2:end) = 0;
                                    theta(2:end) = 0;
                                end
                                anchorPoints = [X; Y; theta];
                                GT_Y = [GT_Y Y-corridorPlannerFrameResampled(indeces,2)];
                                GT_theta = [GT_theta theta];
                                anchorPoints_array = [anchorPoints_array anchorPoints(5:8)-corridorPlannerFrameResampled(indeces,2)];
                                GT_U = [GT_U U];
                            case 4
                                anchorPoints = [X; Y; theta];
                                c0 = refPlannerFrame(:,2)-corridorPlannerFrameResampled(:,2);
                                s = sum((diff(corridorPlannerFrameResampled(:,1)).^2+diff(corridorPlannerFrameResampled(:,2)).^2).^0.5);
                                for i=1:3
                                    clp(i,1) = trapz((c0(indeces(i):indeces(i+1))))/s;
                                end
                                GT_Y = [GT_Y clp];
                                GT_theta = [GT_theta theta];
                                anchorPoints_array = [anchorPoints_array anchorPoints(5:8)-corridorPlannerFrameResampled(indeces,2)];
                                GT_U = [GT_U U];
                        end

                        [trajPlannerFrame, segmentsParam] = trajectory_planner(anchorPoints, indeces, refPlannerFrame,0);
                        trajectoryPreviousPlannerFrame = trajPlannerFrame; %returning previous trajectory in local frame

                        % 7. transforming local trajectories to global one
                        traj_global = ego2GPS(plannerFrame,trajPlannerFrame);
                        ref_global = ego2GPS(plannerFrame,refPlannerFrame);
                        corridor_global = ego2GPS(plannerFrame,corridorPlannerFrameResampled(:,1:2));
                        % 8. Adding trajectories to the globally followed path
                        if (endIndex(1,1)+i > endIdx)
                            traj = [traj; traj_global(:,:)];
                            ref = [ref; ref_global(1:size(traj_global,1),:)];
                            cor = [cor; corridor_global(1:size(traj_global,1),:)];
                            curv = [curv; corridorPlannerFrameResampled(1:size(traj_global,1),4)];
                            orient = [orient; corridorPlannerFrameResampled(1:size(traj_global,1),3)+plannerFrame(3)];
                            validPoints(i:i+size(traj_global,1)-1) = 0;
                            GT_U_array = [GT_U_array zeros(size(GT_U,1), size(traj_global,1))];
                            GT_Y_array = [GT_Y_array zeros(size(GT_Y,1), size(traj_global,1))];
                            break;
                        else
                            traj = [traj; traj_global(counter,:)];
                            ref = [ref; ref_global(counter,:)];
                            cor = [cor; corridor_global(counter,:)];
                            curv = [curv; corridorPlannerFrameResampled(counter,4)];
                            orient = [orient; corridorPlannerFrameResampled(counter,3)+plannerFrame(3)];
                            kappa_array = [kappa_array; mean(GT_U(1:3,end))];
                            validPoints(i) = 1;
                            GT_U_array = [GT_U_array GT_U(:,end)];
                            GT_Y_array = [GT_Y_array GT_Y(:,end)];
                        end
                        if (counter > replanCycle)
                            replan = 2;
                        else
                            replan = 0;
                        end
                    end
                else
                    % some break in the measurement is found, re-init
                    % everything
                    % TODO: check if this if/else is needed, maybe
                    % duplicated
                    replan = 4;
                    initialized = 0;
                    counter = 1; % counter, how far we went in the previous trajectory
                end
            else
                if (i < n)
                    % some break in the measurement is found, re-init
                    % everything
                    replan = 3;
                    initialized = 0;
                    counter = 1; % counter, how far we went in the previous trajectory
                    endIndex =  1000;
                else
                    % we possibly reached the end of the measurement, break the
                    % for loop
                    disp('Measurement end is reached, stopping simulation');
                    break;
                end
            end
        else
            % if there is no replan, the previous trajectory is taken over
            counter = counter + 1;
            if (counter > size(traj_global,1))
                replan = 1;
            else
                if size(traj_global,1)>1
                    traj = [traj; traj_global(counter,:)];
                    ref = [ref; ref_global(counter,:)];
                    cor = [cor; corridor_global(counter,:)];
                    curv = [curv; corridorPlannerFrameResampled(counter,4)];
                    orient = [orient; corridorPlannerFrameResampled(counter,3)+plannerFrame(3)];
                    kappa_array = [kappa_array; GT_U(2,end)];
                    validPoints(i) = 1;
                    GT_U_array = [GT_U_array GT_U(:,end)];
                    GT_Y_array = [GT_Y_array GT_Y(:,end)];
                end

                if (counter > replanCycle)
                    replan = 1;
                else
                    replan = 0;
                end
            end
        end
        %replan_array(i) = replan;
        if (isempty(GT_U))
            %GT_U_array(:,i) = zeros(7,1);
            %GT_Y_array(:,i) = zeros(4,1);
            trajFull(i,:) = [0 0];
            corFull(i,:) = [0 0];
            refFull(i,:) = [0 0];
        else
            %GT_U_array(:,i) = GT_U(:,end);
            %GT_Y_array(:,i) = GT_Y(:,end);
            trajFull = [trajFull; traj(end,:)];
            corFull = [corFull; cor(end,:)];
            refFull = [refFull; ref(end,:)];
        end
    end % END of iterations
    
    %% Postprocessing
    if (isempty(traj))
        corLeft = [];
        corRight = [];
    else
%         valid_indeces = sparePointsFilter(traj);
%         traj = traj(valid_indeces,:);
%         ref = ref(valid_indeces,:);
%         cor = cor(valid_indeces,:);
%         curv = curv(valid_indeces,:);
%         orient = orient(valid_indeces,:);
        % Lane edges
        corLeft = cor + 1.875*ones(size(cor,1),2).*[-sin(orient) cos(orient)];
        corRight = cor + 1.875*ones(size(cor,1),2).*[sin(orient) -cos(orient)];

    end
end

