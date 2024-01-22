function functional_laneWanderingMpc(segments, config)

global snippetsNegOffs i interventionPoint snippetLocalTrajectory snippetLocalVehiclePath snippetLength parametersMpc

RESIMULATE_OPTIMIZED_PARAMETERS = false;
MAXIMUM_NUMBER_OF_SNIPPETS = 6;

temp_folder_path = config.root;
plots_folder_name = 'plots';
set(0,'DefaultFigureVisible','off');

if ~isfolder(fullfile(temp_folder_path, "snippets"))
    mkdir(fullfile(temp_folder_path, "snippets"));
else
    rmdir(fullfile(temp_folder_path, "snippets"),'s');
    mkdir(fullfile(temp_folder_path, "snippets"));
end

optParamFolder = dir(fullfile("C:\git\KDP_Igneczi_new\KDP_Igneczi\publikációk\LaneWandering\data\mpcOptimization\optimizedParametersLeftSuccessful\*.mat"));

for fileID=1:size(segments.segments,2)
    segment = segments.segments(fileID).segment;
    name = segments.segments(fileID).name;
    signalStatus = segments.segments(fileID).signalStatus;
    % writing the metadata for file save at the end
    names{fileID,1} = name(1:end-4);
    
    optParams = [];

    for i=1:length(optParamFolder)
        if (contains(optParamFolder(i).name, name(1:end-8)))
            optParams = load(fullfile(optParamFolder(i).folder, optParamFolder(i).name));
            optParams = optParams.OptimizedParameters;
            break;
        end
    end

    [~, segment_m, indexes] = prepareInputForPlanner(segment);
    curveTypes = calculateCurveType(segment_m, indexes);

    [P0_ma, ~, ~] = functional_calculateStraightLineOffset(segment_m, indexes, 0.111, mean(diff(segment_m(:,indexes.q_T0))));

    offsetError = -segment_m(:,indexes.c0)-P0_ma;
       
    parameters.offsetErrorThd = std(offsetError(curveTypes==1))*0.5;
    parameters.id = name;
    parameters.indexes = indexes;
    
    [snippetsPosOffs, snippetsNegOffs] = functional_laneWanderingSnippetting(offsetError, P0_ma, segment_m, curveTypes, parameters);

    % now rerunning the snippets with proper vehicle init state and calling mpc
    % afterwards to generate the wandering path for both the stray-away and the
    % compensate
    % prediction horizon is estimated based on compensation length
    compensationLength = 0;
    snippetsNegOffs = snippetsPosOffs;

    for i=1:length(snippetsNegOffs)
        compensationLength = compensationLength+(snippetsNegOffs(i).relevantPoints.snippetStop - snippetsNegOffs(i).relevantPoints.snippetStart);
    end

    %% resimulation
    parametersMpc.Np = 45;
    parametersMpc.Nc = 45;

    parametersMpc.mode = 2; % 1 snippet - 1 calibration
   
    % initial parameters for optimization
    x0 = [0; 100; 0; 100];

    if (RESIMULATE_OPTIMIZED_PARAMETERS)
        numberOfSimulatedSnippets = min(MAXIMUM_NUMBER_OF_SNIPPETS,min(size(optParams,1), length(snippetsNegOffs)));
    else
        numberOfSimulatedSnippets = min(MAXIMUM_NUMBER_OF_SNIPPETS, length(snippetsNegOffs));
    end
    
   for i=1:numberOfSimulatedSnippets
       
        %% command line communication
        clc;
        fprintf("mpc resimulation\n");
        fprintf("Simulation/Optimization of driver %d / %d\n", fileID, size(segments.segments,2));
        fprintf("Snippet: %d / %d\n", i, numberOfSimulatedSnippets);

        dataSize = snippetsNegOffs(i).relevantPoints.snippetStop - snippetsNegOffs(i).relevantPoints.snippetStart;
        localOffsetIndex = max(1,snippetsNegOffs(i).relevantPoints.offsetExtremaLocation-snippetsNegOffs(i).relevantPoints.snippetStart);
        
        clear trajectory corridor;
        for j=1:size(snippetsNegOffs(i).snippet,1)
            trajectory(j,1:2) = pos_tf2GPS(...
                snippetsNegOffs(i).snippet(j, snippetsNegOffs(i).indexes.X_abs), ...
                snippetsNegOffs(i).snippet(j, snippetsNegOffs(i).indexes.Y_abs), ...
                snippetsNegOffs(i).snippet(j, snippetsNegOffs(i).indexes.theta_calc), ...
                snippetsNegOffs(i).snippet(j, snippetsNegOffs(i).indexes.c0) + snippetsNegOffs(i).P0_ma(j));
            corridor(j,1:2) = pos_tf2GPS(...
                snippetsNegOffs(i).snippet(j, snippetsNegOffs(i).indexes.X_abs), ...
                snippetsNegOffs(i).snippet(j, snippetsNegOffs(i).indexes.Y_abs), ...
                snippetsNegOffs(i).snippet(j, snippetsNegOffs(i).indexes.theta_calc), ...
                -snippetsNegOffs(i).snippet(j, snippetsNegOffs(i).indexes.c0));
        end

        snippetLocalTrajectory = trajectory - [snippetsNegOffs(i).snippet(1,indexes.X_abs) snippetsNegOffs(i).snippet(1,indexes.Y_abs)];
        T = [cos(snippetsNegOffs(i).snippet(1,indexes.theta_calc)) sin(snippetsNegOffs(i).snippet(1,indexes.theta_calc)); ...
            -sin(snippetsNegOffs(i).snippet(1,indexes.theta_calc)) cos(snippetsNegOffs(i).snippet(1,indexes.theta_calc))];
        snippetLocalTrajectory = snippetLocalTrajectory*T'; % rotational transformation

        snippetLocalVehiclePath = [snippetsNegOffs(i).snippet(:,indexes.X_abs) snippetsNegOffs(i).snippet(:,indexes.Y_abs)] - [snippetsNegOffs(i).snippet(1,indexes.X_abs) snippetsNegOffs(i).snippet(1,indexes.Y_abs)];
        snippetLocalVehiclePath = snippetLocalVehiclePath*T'; % rotational transformation

        snippetCorridor = corridor - [snippetsNegOffs(i).snippet(1,indexes.X_abs) snippetsNegOffs(i).snippet(1,indexes.Y_abs)];
        snippetCorridor = snippetCorridor*T'; % rotational transformation
        
        [snippetsNegOffs(i).snippet, snippetsNegOffs(i).indexes] = addLocalPath(snippetsNegOffs(i).snippet, snippetsNegOffs(i).indexes);

        interventionPoint =  snippetsNegOffs(i).relevantPoints.offsetExtremaLocation - snippetsNegOffs(i).relevantPoints.snippetStart;
        snippetLength = snippetsNegOffs(i).relevantPoints.snippetStop - snippetsNegOffs(i).relevantPoints.snippetStart;

        if (~RESIMULATE_OPTIMIZED_PARAMETERS)
            
            options.OutputFcn = @custom_stop_fun;
            %**************** OPTIMIZATION OF WEIGHTS ******************
            %[x, f] = fmincon(@objectiveFcn,x0,-eye(4),zeros(4,1),[],[],[],[],[]);
            [x, f] = fminsearch(@objectiveFcn, x0);

            fprintf('Optimization of iteration %d has run\n',i);
            fprintf('Final cost is %f, final parameters are [%f %f %f %f]\n', f, x);

            %x0 = x; % returning previous cycle value

            parametersMpc.qDrift = [x(1) 0];
            parametersMpc.rDrift = x(2);
            parametersMpc.qCompensate = [x(3) 0];
            parametersMpc.rCompensate = x(4);

            y = [x; f];

            save(fullfile(temp_folder_path, "snippets",...
                            strcat('OptimizedParameters_', snippetsNegOffs(i).name(1:end-4),'_', num2str(i), '.mat')), "y");
        else
            % ****************** RESIMULATION WITH SPECIFIC WEIGHTS*******
            
            parametersMpc.qDrift = [optParams(i,1) 0];
            parametersMpc.rDrift = optParams(i,2);
            parametersMpc.qCompensate = [optParams(i,3) 0];
            parametersMpc.rCompensate = optParams(i,4);
        end
            
        snippetsNegOffs(i).simulatedPos = ...
            resimulateSnippet(snippetsNegOffs(i).snippet, ...
            snippetsNegOffs(i).indexes, ...
            interventionPoint, ...
            snippetLocalTrajectory, ...
            snippetLocalVehiclePath, ...
            snippetLength, ...
            parametersMpc);

        % calculate the error between the planned and original paths
        ref(:,1) = snippetsNegOffs(i).simulatedPos(:,1);
        ref(:,2) = spline(snippetLocalVehiclePath(:,1), snippetLocalVehiclePath(:,2), ref(:,1));

        if (RESIMULATE_OPTIMIZED_PARAMETERS)
            if (all(optParams(i,:)==0))
                nrms = 'nan';
            else
                e = snippetsNegOffs(i).simulatedPos(:,2) - ref(:,2);
                rms = sum(e.^2)/snippetLength;
                nrms = rms / (max(snippetsNegOffs(i).simulatedPos(:,2))-min(snippetsNegOffs(i).simulatedPos(:,2)));
            end
        else
            e = snippetsNegOffs(i).simulatedPos(:,2) - ref(:,2);
            rms = sum(e.^2)/snippetLength;
            nrms = rms / (max(snippetsNegOffs(i).simulatedPos(:,2))-min(snippetsNegOffs(i).simulatedPos(:,2)));
        end


        data(i).snippetLocalTrajectory = snippetLocalTrajectory;
        data(i).snippetLocalVehiclePath = snippetLocalVehiclePath;
        data(i).snippetsNegOffs = snippetsNegOffs(i);
        data(i).e = nrms;
        clear ref;

    % ******************* PLOTTING **********************************
        if (size(snippetsNegOffs(i).simulatedPos,1) > 5)
    
            f = figure();

            set(f,'defaulttextInterpreter','latex') ;
            set(f, 'defaultAxesTickLabelInterpreter','latex');  
            set(f, 'defaultLegendInterpreter','latex');
    
            subplot(2,1,1);
            % path output
            plot(snippetsNegOffs(i).simulatedPos(:,1), snippetsNegOffs(i).simulatedPos(:,2), 'r--', 'LineWidth', 2);
            hold on; grid on;
            plot(snippetLocalTrajectory(1:dataSize,1), snippetLocalTrajectory(1:dataSize,2), 'k');
            plot(snippetLocalVehiclePath(1:dataSize,1), snippetLocalVehiclePath(1:dataSize,2), 'b');
            xline(snippetLocalTrajectory(localOffsetIndex), ...
                'color', 'k', 'LineWidth', 2);
    
    
            lgd = legend('Simulation', 'Path', 'Human', 'Intervention Point');
            lgd.Orientation = 'vertical';
            lgd.Location = 'best';
            set(gca,'FontSize', 12);
            xlabel('x(m)', 'FontSize', 12); ylabel('y(m)', 'FontSize', 12);            
            title('\textbf{Path selection during drifting and compensation}', 'FontSize',14);
            xlim([min(snippetLocalTrajectory(1:dataSize,1)), max(snippetLocalTrajectory(1:dataSize,1))]);
    
            subplot(2,1,2);
            % offset output
            % simulated output resampled at X of reference
            pathY = spline(snippetsNegOffs(i).simulatedPos(:,1), snippetsNegOffs(i).simulatedPos(:,2), snippetLocalTrajectory(1:dataSize,1));
            simulatedOffset = pathY - snippetLocalTrajectory(1:dataSize,2);
            plot(snippetLocalTrajectory(1:dataSize,1), snippetsNegOffs(i).offsetError(1:dataSize), 'b');
            hold on; grid on;
                
            plot(snippetLocalTrajectory(1:dataSize,1), simulatedOffset, 'r--', 'LineWidth',2);
            ylim([-0.5, 0.5]);
            xlim([min(snippetLocalTrajectory(1:dataSize,1)), max(snippetLocalTrajectory(1:dataSize,1))]);
            xline(snippetLocalTrajectory(localOffsetIndex), ...
                'color', 'k', 'LineWidth', 2);
    
            set(gca,'FontSize', 12);
            xlabel('x(m)', 'FontSize', 12);
            ylabel('offset error(m)', 'FontSize', 12);
            lgd = legend('Human', 'Simulated');
            lgd.Orientation = 'horizontal';
            lgd.Location = 'south';
            title('\textbf{Offset error}', 'FontSize',14);
    
            savefig(f, fullfile(temp_folder_path, "snippets",...
                            strcat('Resimulated_snippets', snippetsNegOffs(i).name(1:end-4),'_', num2str(i), '.fig')));
            saveas(f, fullfile(temp_folder_path, "snippets",...
                            strcat('Resimulated_snippets', snippetsNegOffs(i).name(1:end-4),'_', num2str(i), '.png')));
            close(f);
        end
   end
   if (numberOfSimulatedSnippets > 0)
        save(fullfile(temp_folder_path, "snippets",...
                             strcat('DataMpc', name(1:end-4), '.mat')), 'data');
   end
   clear data;
end
  
end   

function stop = custom_stop_fun(~, optimValues, ~)
if optimValues.iteration > 200
    stop = true;
elseif optimValues.iteration > 40
        if optimValues.fval <=3e-3
            stop = true;
        else
            stop = false;
        end
elseif optimValues.iteration > 3 && optimValues.fval > 2
    stop = true;
elseif optimValues.fval <=1e-3
    stop = true;
else
    stop = false;
end
clc;
fprintf("Optimization cycle is %d\n", optimValues.iteration);
fprintf("Cost is %f\n", optimValues.fval);
end

function f = objectiveFcn(x0)
global snippetsNegOffs i interventionPoint snippetLocalTrajectory snippetLocalVehiclePath snippetLength parametersMpc

    parametersMpc.qDrift = [x0(1) 0];
    parametersMpc.rDrift = x0(2);
    parametersMpc.qCompensate = [x0(3) 0];
    parametersMpc.rCompensate = x0(4);

    [simulatedPos] = resimulateSnippet(snippetsNegOffs(i).snippet, snippetsNegOffs(i).indexes, ...
                interventionPoint, ...
                snippetLocalTrajectory, ...
                snippetLocalVehiclePath, ...
                snippetLength, ...
                parametersMpc);

    dataSize = size(snippetLocalVehiclePath,1);

    ref(:,1) = simulatedPos(:,1);
    ref(:,2) = spline(snippetLocalVehiclePath(:,1), snippetLocalVehiclePath(:,2), ref(:,1));

%     p = plot(ref(:,1), ref(:,2), simulatedPos(:,1), simulatedPos(:,2), snippetLocalTrajectory(:,1), snippetLocalTrajectory(:,2), snippetLocalVehiclePath(interventionPoint,1), snippetLocalVehiclePath(interventionPoint,2));
%     set(p(3), 'color', 'k', 'LineStyle', '--');
%     set(p(4), 'color', 'k', 'Marker', 'o', 'LineWidth', 3);
%     shg;

    f = sum((simulatedPos(:,2) - ref(:,2)).^2);
    
    f = f/dataSize + 1000*(-1)*min(min(0, x0));
end

function simulatedPos = resimulateSnippet(snippet, indexes, interventionPoint, snippetLocalTrajectory, snippetLocalVehiclePath, snippetEnd, parameters)
    
    localOffsetIndex = max(1,interventionPoint);

    vehicleState = initVehicleState(snippet, indexes, 1, snippetLocalVehiclePath);
    dT = mean(diff(snippet(:,indexes.q_T0)));
    
    x_k1 =  [0; 0; vehicleState.vx; 0; 0];
    j = 1;
    run = true;
    interventionPoint = localOffsetIndex;
    while (run)
        simulatedPos(j,1:2) = [vehicleState.X vehicleState.Y];

        dX = vehicleState.vx*0.02;
        x = linspace(0,dX*(parameters.Np-1), parameters.Np);

        % cutting the global path to start from the nearest point
        % find nearest point
        idx = find(snippetLocalTrajectory(:,1) >= vehicleState.X,1);
        if (isempty(idx))
            break;
        else
            if ((snippetLocalTrajectory(end,1) - snippetLocalTrajectory(idx(1,1),1)) < x(end))
                break;
            end
        end
        if (vehicleState.X > snippetLocalTrajectory(snippetEnd, 1))
            % the end of the actual compensation phase is reached
            break;
        end
        if (j>1000)
            break;
        end

        % resampling path for prediction points
        x = x+snippetLocalTrajectory(idx(1,1),1);
        y = spline(snippetLocalTrajectory(:,1), snippetLocalTrajectory(:,2), x);

        pathLocal = [x' y'];
        clear x y
        pathOrientation = diff(pathLocal(:,2))./diff(pathLocal(:,1));
        pathOrientation = [pathOrientation; pathOrientation(end)];
        pathOrientation = movmean(pathOrientation, 25);
        pathOrientation = atan(pathOrientation);

        x = [vehicleState.X; vehicleState.Y; vehicleState.vx; vehicleState.vy; vehicleState.theta];
        xa = [x-x_k1; vehicleState.Y; vehicleState.theta];        
        
        % adaptive parametrization
        if (parameters.mode == 1)
            q = parameters.q;
            rw = parameters.rw; 
        elseif (parameters.mode == 2)
            if (vehicleState.X < snippetLocalTrajectory(localOffsetIndex,1))
                % before interaction point
                q = parameters.qDrift;
                rw = parameters.rDrift; 
            else
                % left drift
                q = parameters.qCompensate;
                rw = parameters.rCompensate;
            end
        end
        ayTar = mpc(pathLocal, pathOrientation, parameters.Np,  parameters.Nc, rw, q, dT,  x_k1, xa, vehicleState.ay_v, vehicleState.theta, vehicleState.vx_v, 2.5, vehicleState.ay_v);
        vehicleState.yawRate = ayTar/vehicleState.vx;
        vehicleState = vehicleModel(vehicleState, 0.02);

        x_k1 = x;
        j=j+1;
    end

end

function [P0_ma, P0, P0_fd, P0_ocNormal, P0_ocTruck, error_fd, error_ocNormal, error_ocTruck, errorValues_fd, errorValues_ocNormal, errorValues_ocTruck, data] = calculateStraightLineOffset(segment_m, indexes, curveTypes, fc)
    error_ocNormal = []; error_ocTruck = []; errorValues_fd= []; errorValues_ocNormal=[];errorValues_ocTruck=[];
    oncomingTrafficNormal = zeros(size(segment_m,1),1);
    oncomingTraffic = zeros(size(segment_m,1),1);
    oncomingTrafficTruck = zeros(size(segment_m,1),1);

    oncomingTraffic(segment_m(:,indexes.oncomingTraffic) > 0) = 1;
    oncomingTrafficNormal(segment_m(:,indexes.oncomingTraffic) == 2 | segment_m(:,indexes.oncomingTraffic) == 22) = 1;
    oncomingTrafficTruck(~oncomingTrafficNormal & oncomingTraffic) = 1;

    data.oncomingTraffic = oncomingTraffic;
    data.oncomingTrafficNormal = oncomingTrafficNormal;
    data.oncomingTrafficTruck = oncomingTrafficTruck;

    oncomingTrafficNormal(curveTypes~=1) = 0;
    oncomingTrafficTruck(curveTypes~=1) = 0;

    oncomingTraffic = morphologyClose(oncomingTraffic, 25);
    oncomingTrafficNormal = morphologyClose(oncomingTrafficNormal, 25);    
    oncomingTrafficTruck = morphologyClose(oncomingTrafficTruck, 25);  

    % a new indicator is introduced, which represents two things:
    % 1) the time since the last vehicle passed us
    % 2) the time until the next vehicle is in sight
    % this number is characterized by a number between 0 and 1. When a
    % vehicle left us, with a certain gradient its psychological effect is
    % downgraded. If there is a suspected vehicle in sight, it raises the
    % awareness of the driver, with a certain gradient counting up.
    % literature: 2-3 s of short term memory.
    risingEdges = find(diff(oncomingTraffic)>0);
    fallingEdges = find(diff(oncomingTraffic)<0);
    oncomingTrafficScaled = oncomingTraffic;
    N = 100;
    risingSections = linspace(0,1,N)';
    fallingSections = risingSections(end:-1:1);

    oncomingTrafficScaledWithRisingEdges = oncomingTrafficScaled;
    for i=1:length(risingEdges)
        oncomingTrafficScaledWithRisingEdges(risingEdges(i)-(N-1):risingEdges(i)) = risingSections;
    end
    oncomingTrafficScaledWithFallingEdges = oncomingTrafficScaled;
    for i=1:length(fallingEdges)
        oncomingTrafficScaledWithFallingEdges(fallingEdges(i):fallingEdges(i)+N-1) = fallingSections;
    end

    oncomingTrafficScaledWithFallingEdges = oncomingTrafficScaledWithFallingEdges(1:size(oncomingTraffic,1),1);

    data.oncomingTrafficScaled = max(oncomingTrafficScaledWithRisingEdges, oncomingTrafficScaledWithFallingEdges);

    risingEdges = find(diff(oncomingTrafficTruck)>0);
    fallingEdges = find(diff(oncomingTrafficTruck)<0);
    oncomingTrafficTruckScaled = oncomingTrafficTruck;
    N = 100;
    risingSections = linspace(0,1,N)';
    fallingSections = risingSections(end:-1:1);

    oncomingTrafficScaledWithRisingEdges = oncomingTrafficTruckScaled;
    for i=1:length(risingEdges)
        oncomingTrafficScaledWithRisingEdges(risingEdges(i)-(N-1):risingEdges(i)) = risingSections;
    end
    oncomingTrafficScaledWithFallingEdges = oncomingTrafficTruckScaled;
    for i=1:length(fallingEdges)
        oncomingTrafficScaledWithFallingEdges(fallingEdges(i):fallingEdges(i)+N-1) = fallingSections;
    end

    oncomingTrafficScaledWithFallingEdges = oncomingTrafficScaledWithFallingEdges(1:size(oncomingTraffic,1),1);

    oncomingTrafficTruck = max(oncomingTrafficScaledWithRisingEdges, oncomingTrafficScaledWithFallingEdges);

    risingEdges = find(diff(oncomingTrafficNormal)>0);
    fallingEdges = find(diff(oncomingTrafficNormal)<0);
    oncomingTrafficNormalScaled = oncomingTrafficNormal;
    N = 100;
    risingSections = linspace(0,1,N)';
    fallingSections = risingSections(end:-1:1);

    oncomingTrafficScaledWithRisingEdges = oncomingTrafficNormalScaled;
    for i=1:length(risingEdges)
        oncomingTrafficScaledWithRisingEdges(risingEdges(i)-(N-1):risingEdges(i)) = risingSections;
    end
    oncomingTrafficScaledWithFallingEdges = oncomingTrafficNormalScaled;
    for i=1:length(fallingEdges)
        oncomingTrafficScaledWithFallingEdges(fallingEdges(i):fallingEdges(i)+N-1) = fallingSections;
    end

    oncomingTrafficScaledWithFallingEdges = oncomingTrafficScaledWithFallingEdges(1:size(oncomingTraffic,1),1);

    oncomingTrafficNormal = max(oncomingTrafficScaledWithRisingEdges, oncomingTrafficScaledWithFallingEdges);

    freeDriving = zeros(size(segment_m,1),1);
    freeDriving(~oncomingTraffic) = 1;
    freeDriving(curveTypes~=1) = 0;

    P0_ocNormal = zeros(size(segment_m,1),1);
    P0_ocTruck= zeros(size(segment_m,1),1);
    P0_fd = zeros(size(segment_m,1),1);

    P0 = zeros(size(segment_m,1),1);
    P0_ma = movmean(-segment_m(:, indexes.c0), window);
    i = 1;
    while(true)
        if (i==1)
            if (oncomingTrafficNormal(i) == 0)
                % start with freedriving
                oncomingEnd = 1;
                oncomingStart = find(oncomingTrafficNormal(oncomingEnd:end)==1,1) + oncomingEnd;
                if (isempty(oncomingStart))
                    break;
                end
                oncomingEnd = find(oncomingTrafficNormal(oncomingStart:end)==0,1)-1+oncomingStart;
                if (isempty(oncomingEnd))
                    break;
                end
                P0_oc(oncomingStart:oncomingEnd) = mean(-segment_m(oncomingStart:oncomingEnd,indexes.c0));
            else
                % start with oncoming
                oncomingStart = 1;
                oncomingEnd = find(oncomingTrafficNormal(oncomingStart:end)==0,1)-1+oncomingStart;
                if (isempty(oncomingEnd))
                    break;
                end
                P0_oc(oncomingStart:oncomingEnd) = mean(-segment_m(oncomingStart:oncomingEnd,indexes.c0));
            end
        else
            oncomingStart = find(oncomingTrafficNormal(oncomingEnd:end)==1,1)-1 + oncomingEnd;
            if (isempty(oncomingStart))
                break;
            end
                
            oncomingEnd = find(oncomingTrafficNormal(oncomingStart:end)==0,1)-1+oncomingStart;
            if (isempty(oncomingEnd))
                    break;
            end
            P0_ocNormal(oncomingStart:oncomingEnd) = mean(-segment_m(oncomingStart:oncomingEnd,indexes.c0));
        end
        error_ocNormal{i} = P0_ocNormal(oncomingStart:oncomingEnd) + segment_m(oncomingStart:oncomingEnd,indexes.c0);
        [hist_ocNormal(i).N, hist_ocNormal(i).edges] = histcounts(error_ocNormal{i});
        errorValues_ocNormal(i) = P0_ocNormal(oncomingStart);
        i = i+1;
    end

    i = 1;
    while(true)
        if (i==1)
            if (oncomingTrafficTruck(i) == 0)
                % start with freedriving
                oncomingEnd = 1;
                oncomingStart = find(oncomingTrafficTruck(oncomingEnd:end)==1,1) + oncomingEnd;
                if (isempty(oncomingStart))
                    break;
                end
                oncomingEnd = find(oncomingTrafficTruck(oncomingStart:end)==0,1)-1+oncomingStart;
                if (isempty(oncomingEnd))
                    break;
                end
                P0_oc(oncomingStart:oncomingEnd) = mean(-segment_m(oncomingStart:oncomingEnd,indexes.c0));
            else
                % start with oncoming
                oncomingStart = 1;
                oncomingEnd = find(oncomingTrafficTruck(oncomingStart:end)==0,1)-1+oncomingStart;
                if (isempty(oncomingEnd))
                    break;
                end
                P0_oc(oncomingStart:oncomingEnd) = mean(-segment_m(oncomingStart:oncomingEnd,indexes.c0));
            end
        else
            oncomingStart = find(oncomingTrafficTruck(oncomingEnd:end)==1,1)-1 + oncomingEnd;
            if (isempty(oncomingStart))
                break;
            end
                
            oncomingEnd = find(oncomingTrafficTruck(oncomingStart:end)==0,1)-1+oncomingStart;
            if (isempty(oncomingEnd))
                    break;
            end
            P0_ocTruck(oncomingStart:oncomingEnd) = mean(-segment_m(oncomingStart:oncomingEnd,indexes.c0));
        end
        error_ocTruck{i} = P0_ocTruck(oncomingStart:oncomingEnd) + segment_m(oncomingStart:oncomingEnd,indexes.c0);
        [hist_ocTruck(i).N, hist_ocTruck(i).edges] = histcounts(error_ocTruck{i});
        errorValues_ocTruck(i) = P0_ocTruck(oncomingStart);
        i = i+1;
    end

    i = 1;
    while(true)
        if (i==1)
            if (freeDriving(i) == 0)
                % start with freedriving
                
                freeDrivingEnd = 1;
                freeDrivingStart = find(freeDriving(freeDrivingEnd:end)==1,1) + freeDrivingEnd;
                if (isempty(freeDrivingStart))
                    break;
                end
                freeDrivingEnd = find(freeDriving(freeDrivingStart:end)==0,1)-1+freeDrivingStart;
                if (isempty(freeDrivingEnd))
                    break;
                end
                P0_fd(freeDrivingStart:freeDrivingEnd) = mean(-segment_m(freeDrivingStart:freeDrivingEnd,indexes.c0));
            else
                % start with oncoming
                freeDrivingStart = 1;
                freeDrivingEnd = find(freeDriving(freeDrivingStart:end)==0,1)-1+freeDrivingStart;
                if (isempty(freeDrivingEnd))
                    break;
                end
                P0_fd(freeDrivingStart:freeDrivingEnd) = mean(-segment_m(freeDrivingStart:freeDrivingEnd,indexes.c0));
            end
        else
            freeDrivingStart = find(freeDriving(freeDrivingEnd:end)==1,1)-1 + freeDrivingEnd;
            if (isempty(freeDrivingStart))
                break;
            end
                
            freeDrivingEnd = find(freeDriving(freeDrivingStart:end)==0,1)-1+freeDrivingStart;
            if (isempty(freeDrivingEnd))
                    break;
            end
            P0_fd(freeDrivingStart:freeDrivingEnd) = mean(-segment_m(freeDrivingStart:freeDrivingEnd,indexes.c0));
        end
        % calculating the histogram data - error between c0 and averaged
        error_fd{i} = P0_fd(freeDrivingStart:freeDrivingEnd) + segment_m(freeDrivingStart:freeDrivingEnd,indexes.c0);
        [hist_fd(i).N, hist_fd(i).edges] = histcounts(error_fd{i});
        errorValues_fd(i) = P0_fd(freeDrivingStart);
        i = i+1;
    end

    P0 = P0_fd+P0_ocTruck+P0_ocNormal;
end

function [segment_m, indexes] = addLocalPath(segment_m, indexes)


theta0 = segment_m(1, indexes.theta_calc);
T = [cos(theta0) sin(theta0); -sin(theta0) cos(theta0)];

refLocal = [segment_m(:,indexes.X_abs) segment_m(:,indexes.Y_abs)] - [segment_m(1,indexes.X_abs) segment_m(1,indexes.Y_abs)];
refLocal = refLocal*T';
segment_m(:,end+1:end+2) = refLocal;

indexes.X = size(segment_m,2)-1;
indexes.Y = size(segment_m,2);

end

function vehicleState = initVehicleState(segment_m, indexes, index, snippetLocalVehiclePath)
    theta = atan(diff(snippetLocalVehiclePath(:,2))./diff(snippetLocalVehiclePath(:,1)));
    %theta = movmean(theta,5);
    theta0 = mean(theta(1:min(size(theta,1),20)));
    p_wheelBase = 2.7;
    vehicleState.X = 0; %segment_m(index, indexes.X_abs);
    vehicleState.Y = 0; %segment_m(index, indexes.Y_abs);
    vehicleState.theta = theta0; %segment_m(index, indexes.theta_calc);
    yawRate = movmean(segment_m(:, indexes.yawRate), 50);
    vehicleState.yawRate = yawRate(index);
    vehicleState.vx_v = segment_m(index, indexes.velocityX)*cos(theta0);
    vehicleState.vy_v = segment_m(index, indexes.velocityX)*sin(theta0);
    vehicleState.ax_v = 0;
    vehicleState.ay_v = vehicleState.vx_v*vehicleState.yawRate;
    vehicleState.t0 = segment_m(index, indexes.q_T0);    
    vehicleState.steeringAngle = 0;

    % initial global states
    vehicleState.vx = vehicleState.vx_v*cos(vehicleState.theta);
    vehicleState.vy = vehicleState.vx_v*sin(vehicleState.theta);
end

function vehicleState = vehicleModel(vehicleState, dT)
    % kinematic bicycle model
    vehicleState.theta = vehicleState.theta + dT*vehicleState.yawRate;
    % displacement
    dX = vehicleState.vx_v*dT*cos(vehicleState.theta);
    dY = vehicleState.vx_v*dT*sin(vehicleState.theta);

    % new vehicle position
    vehicleState.X = vehicleState.X+dX;
    vehicleState.Y = vehicleState.Y+dY;

    % new global states
    vehicleState.vx = vehicleState.vx_v*cos(vehicleState.theta);
    vehicleState.vy = vehicleState.vx_v*sin(vehicleState.theta);
    vehicleState.ax_v = 0;
    vehicleState.ay_v = vehicleState.vx_v*vehicleState.yawRate;    
end

function [u, Y] = mpc(pathLocal, pathOrientation, Np, Nc, rw, q, Ts,  x_k1, xa, aeta0, theta0, v0, ay_max, u_k1)
%% INTRODUCTION
% This function is the core MPC algorithm part.
% created by Gergo Igneczi @ Vehicle Research Center of Szechenyi Istvan
% University

if (x_k1(3) < 3)
    % speed is low, MPC will not work
    u = 0;
else
    dim = 2; % number of outputs
    % producing prediction matrices
    pred_matrix = zeros(min(dim*Np,1000),1);

    x_points = pathLocal(:,1);
    y_points = pathLocal(:,2);

    for i = 1:Np
        %pred_matrix(dim*i-(dim-1),1) = x_points(i);
        pred_matrix(dim*i-(dim-1),1) = y_points(i);
        pred_matrix(dim*i-(dim-2),1) = pathOrientation(i);
    end
    Rs_rk = pred_matrix;

    I = eye(min(Nc,1000));
    %% determining state constraints    

    u_max = ay_max;
    u_min =  ay_max;

    %% initializing helper matrices
    M = zeros(min(1000,2*size(u_k1,1)*Nc),min(2*Nc,1000));
    for (i=1:2*size(u_k1,1)*Nc)
        if (i<=(Nc))
            k=1;
            while ((2*k-1)/2 <= i)
                M (i,2*k-1) = 1;
                k = k + 1;
            end
        elseif (i <=(Nc*2))
            k=1;
            while ((2*k-1)/2 <= i-Nc)
                M (i,2*k-1) = -1;
                k = k + 1;
            end
        elseif (i<=(Nc*3))
            k=1;
            while (2*k/2 <= i-size(x_k1,1)*Nc/2)
                M (i,2*k) = 1;
                k = k+1;
            end
        else
            k=1;
            while (2*k/2 <= i-size(x_k1,1)*Nc/2-Nc)
                M (i,2*k) = -1;
                k = k+1;
            end
        end
    end


    %% updating state matrices
    %Ad = [1 0 Ts 0 0; 0 1 0 Ts 0; 0 0 1 0 -Ts/L*x_k1(3)^2*tan(delta0); 0 0 0 1 0; 0 0 0 0 1];
    %Bd = [0 0 0 Ts/L*x_k1(3)^2*1/(cos(delta0))^2 Ts/L*x_k1(3)*1/(cos(delta0))^2]';
    Ad = [1 0 Ts 0 0; 0 1 0 Ts 0; 0 0 1 0 -aeta0*Ts*cos(theta0); 0 0 0 1 -Ts*aeta0*sin(theta0);0 0 0 0 1];
    Bd = [0 0 -Ts*sin(theta0) Ts*cos(theta0) Ts*aeta0/v0]';
    Cd = [0 1 0 0 0; 0 0 0 0 1];
    
    n = size(Ad,1); m = size(Cd,1); k = 1; %size(Bd,2);
    
    %% augmented model
    A = [Ad zeros(n,m); Cd*Ad eye(m,m)];
    B = [Bd; Cd*Bd];
    C = [zeros(m,n) eye(m,m)];
    %% matrix generation
    F = zeros(Np*m,m+n);
    for i=1:Np
        if (m>1)
            F(m*i-(m-1):m*i,1:(m+n))=C*A^i;
        else
            F(i,1:(m+n))=C*A^i;
        end
    end
    S = zeros(m*Np,Nc);
    for i=1:Np
        for j=1:Nc
            if(j>i)
                S(m*i-(m-1):m*i,j)=zeros(m,k);
            else
                S(m*i-(m-1):m*i,j)=C*A^((i-1)-(j-1))*B;
            end
        end
    end

    %% control calculation
    % Unconstrained results
    dU = zeros(min(1000,k*Nc),1);
    %dU = inv(S'*S+rw*I)*(S'*Rs_rk-S'*F*xa);
    R = rw*I;
    Q = zeros(m*Np, m*Np);
    for i=1:Np
        Q(i*numel(q)-(numel(q)-1):i*numel(q),i*numel(q)-(numel(q)-1):i*numel(q)) = diag(q);
    end
    dU = inv(S'*Q*S+R)*S'*Q*(Rs_rk-F*xa);
 
    % Constrained results
    gamma = zeros(min(1000,2*size(u_k1,1)*Nc),1);

    for i=1:2*size(u_k1,1)*Nc
        if (i<=Nc)
            gamma(i,1) = u_max(1)-u_k1(1);
        elseif (i<=2*Nc)
            gamma(i,1) = u_min(1)+u_k1(1);
        end
    end

    if (all(M(:,1:2:end)*dU<=gamma))
        %do nothing
    else
        %Solving Hildreth's QP problem
        E = (S'*S+rw*I)*2;
        F_ = -2*S'*(Rs_rk-F*xa);
        H = E;
        f = F_;
        A_cons = M(:,1:2:end);
        b = gamma;
        eta = x_k1;
        [n1,m1]=size(A_cons);
        eta=-H\f;
        kk=0;
        for i=1:n1
            if (A_cons(i,:)*eta>b(i)) 
                kk=kk+1;
            else
                kk=kk+0;
            end
        end
        if (kk==0) 
            % do nothing 
        else
            P=A_cons*(H\A_cons');
            d=(A_cons*(H\f)+b);
            [n,m]=size(d);
            x_ini=zeros(n,m);
            lambda=x_ini;
            al=10;
            for km=1:38
                %find the elements in the solution vector one by one
                % km could be larger if the Lagranger multiplier has a slow
                % convergence rate.
                lambda_p=lambda;
                for i=1:n
                    w= P(i,:)*lambda-P(i,i)*lambda(i,1);
                    w=w+d(i,1);
                    la=-w/P(i,i);
                    lambda(i,:)=max(0,la);
                end
                al=(lambda-lambda_p)'*(lambda-lambda_p);
                if (al<10e-8)
                    break; 
                end
            end
            dU=-H\f -H\A_cons'*lambda;
        end
    end

    u = u_k1+dU(1);
    Y = F*xa+S*dU;
end
end
