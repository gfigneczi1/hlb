function evaluator_laneWandering(segments,config)

% This evaluator is responsible for generating various plots to analyse
% lane wandering performance of different drivers.
% The following KPI-s are proposed to be used:
% - straight line:
% -- compensation points (where the steering torque is applied to compensate stray effect)
% -- distribution of lane offset and CLP
% -- oncoming traffic? follow traffic?

% Segmentor: segmentor_driverModel.

% @copyright (C) 2022 Robert Bosch GmbH.
% The reproduction, distribution and utilization of this file as
% well as the communication of its contents to others without express
% authorization is prohibited. Offenders will be held liable for the
% payment of damages. All rights reserved in the event of the grant
% of a patent, utility model or design.
% @version 1.0
% @date 2023-09-04

% step 1: interaction point model relations
%interactionModelling (segments, config);
%offsetPlotting(segments, config);

% step 2: snippeting
% cutting the snipets where the extremum is detected until the point, where
% vehicle settles around the average offset
% 2.1 snippet statistics - either ttcl or offset based
params.edgesX_Left = linspace(0,0.8,10);
params.edgesX_Right = linspace(-0.8,0,10);
params.edgesX_Left = linspace(0,10,10);
params.edgesX_Right = linspace(0,10,10);
params.xlabelstr = "offset error(m)";

%functional_laneWanderingDistribution(segments, config, params);

%snippets = cutSnippets(segments, config, params);

%snippetStatisticsOffsetBased (segments, config, params);

% step 4: Running MPC for snippets
windowChecks(segments, config, params);
%snippetModelByMpc(segments, config, params);

end

function [data] = cutSnippets(segments, config, params)

temp_folder_path = config.root;
plots_folder_name = 'plots';
set(0,'DefaultFigureVisible','off');

if ~isfolder(fullfile(temp_folder_path, "snippets"))
    mkdir(fullfile(temp_folder_path, "snippets"));
else
    rmdir(fullfile(temp_folder_path, "snippets"),'s');
    mkdir(fullfile(temp_folder_path, "snippets"));
end

for fileID=1:size(segments.segments,2)
    segment = segments.segments(fileID).segment;
    name = segments.segments(fileID).name;
    signalStatus = segments.segments(fileID).signalStatus;
    % writing the metadata for file save at the end
    names{fileID,1} = name(1:end-4);

    [~, segment_m, indexes] = prepareInputForPlanner(segment);
    curveTypes = calculateCurveType(segment_m, indexes, name);
    window = linspace(50, 500, 11); % with step size of 50 samples = 2.5s

    [P0_ma, P0, P0_fd, P0_oc, error_fd, error_oc, errorValues_fd, errorValues_oc] = calculateStraightLineOffset(segment_m, indexes, curveTypes, 180);

    offsetError = -segment_m(:,indexes.c0)-P0_ma;
    thd = std(offsetError(curveTypes==1))*0.5;
    interactionPoints = abs(offsetError) > thd;
    interactionPoints(curveTypes ~=1) = false;
    interactionPoints = morphologyOpen(interactionPoints, 10);
    interactionPoints = morphologyClose(interactionPoints, 10);
    extrema = findExtremumPoints(offsetError, interactionPoints);

    %% cut snippets
    j = 1;
    k = 1;
    kpisNegOffs = [];
    kpisPosOffs = [];
    %clear snippetsNegOffs snippetsPosOffs
    snippetsNegOffs = [];
    snippetsPosOffs = [];

    for i=1:size(extrema,1)
        snippetStart = extrema(i,1);
        % snippet start is modified to get the last point, where the thd
        % was exceeded before reaching the extremum
        offsetsBackwards = offsetError(snippetStart:-1:1);
        if (extrema(i,2) > 0)
            % this is a positive offset, searching for first point
            % backwards with negative sign
            thresholdExceed = find(offsetsBackwards <= 0, 1);
        else
            thresholdExceed = find(offsetsBackwards >= 0, 1);
        end

        snippetStop = snippetStart+find(abs(offsetError(snippetStart:end)) < 0.02,1); % the first point where c0 reaches its average after extremum is detected
        if (~isempty(thresholdExceed))
            snippetStart = snippetStart-thresholdExceed;
        end
        if (isempty(snippetStop))
           break;
        end
        compensationAborted = false;
        if (i<size(extrema,1))
            if (snippetStop > extrema(i+1,1))
                % this is the situation, where something happened before
                % getting back to the settle offset, therefore we simply
                % neglect this compensation
                compensationAborted = true;
            end
        end
        if (compensationAborted)
        else
            snippetStop = min(snippetStop, size(segment_m,1));
            snippetStart = max(snippetStart,1); % adding some init phase
            snippetStopExtended = min(snippetStop + (snippetStop-snippetStart)*2, size(segment_m,1));
            if ((snippetStop-snippetStart) <= 180)
                if (extrema(i,2) > 0)
                    % the offset extremum is above the straight line offset
                    % of the driver
                    snippetsPosOffs(j).name = name;
                    snippetsPosOffs(j).snippet = segment_m(snippetStart:snippetStopExtended, :);
                    snippetsPosOffs(j).indexes = indexes;
                    snippetsPosOffs(j).offsetError = offsetError(snippetStart:snippetStopExtended);
                    snippetsPosOffs(j).P0_ma = P0_ma(snippetStart:snippetStopExtended);
                    snippetsPosOffs(j).relevantPoints.snippetStart = snippetStart;
                    snippetsPosOffs(j).relevantPoints.snippetStop = snippetStop;
                    snippetsPosOffs(j).relevantPoints.snippetStopExtended = snippetStopExtended;
                    snippetsPosOffs(j).relevantPoints.offsetExtremaLocation = extrema(i,1);
                    % for positive offset the following quantities are
                    % calculated&saved
                    % - minimum absolute relative lateral acceleration
                    % (right side acceleration)
                    % - minimum absolute relative lateral jerk
                    % (right side jerk)
                    % - minimum absolute curvature change (right side
                    % curvature gradient)
                    % - maximum lateral offset to the left side
                    kpisPosOffs(j, 1:5) = [min(segment_m(snippetStart:snippetStop, indexes.ayRel)) ...
                        min(segment_m(snippetStart:snippetStop, indexes.jyRel)) ...
                        min(segment_m(snippetStart:snippetStop, indexes.vehicleCurvatureChange)) ...
                        min(segment_m(snippetStart:snippetStop, indexes.TTCL)) ...
                        max(offsetError(snippetStart:snippetStop))];
                     j= j+1;
                else
                    snippetsNegOffs(k).name = name;
                    snippetsNegOffs(k).snippet = segment_m(snippetStart:snippetStopExtended, :);
                    snippetsNegOffs(k).indexes = indexes;
                    snippetsNegOffs(k).offsetError = offsetError(snippetStart:snippetStopExtended);
                    snippetsNegOffs(k).P0_ma = P0_ma(snippetStart:snippetStopExtended);
                    snippetsNegOffs(k).relevantPoints.snippetStart = snippetStart;
                    snippetsNegOffs(k).relevantPoints.snippetStop = snippetStop;
                    snippetsNegOffs(k).relevantPoints.snippetStopExtended = snippetStopExtended;
                    snippetsNegOffs(k).relevantPoints.offsetExtremaLocation = extrema(i,1);
                    kpisNegOffs(k, 1:5) = [max(segment_m(snippetStart:snippetStop, indexes.ayRel)) ...
                        max(segment_m(snippetStart:snippetStop, indexes.jyRel)) ...
                        max(segment_m(snippetStart:snippetStop, indexes.vehicleCurvatureChange)) ...
                        min(segment_m(snippetStart:snippetStop, indexes.TTCL)) ...
                        min(offsetError(snippetStart:snippetStop))];
                    k = k+1;
                end
            end
        end
    end 

    % now rerunning the snippets with proper vehicle init state and calling mpc
    % afterwards to generate the wandering path for both the stray-away and the
    % compensate
    % prediction horizon is estimated based on compensation length
    compensationLength = 0;

    for i=1:length(snippetsNegOffs)
        compensationLength = compensationLength+(snippetsNegOffs(i).relevantPoints.snippetStop - snippetsNegOffs(i).relevantPoints.snippetStart);
    end

    data(fileID).snippetsNegOffs = snippetsNegOffs;
    data(fileID).snippetsPosOffs = snippetsPosOffs;
    data(fileID).kpisNegOffs = kpisNegOffs;
    data(fileID).kpisPosOffs = kpisPosOffs;
    data(fileID).name = name;
end

end

function windowChecks(segments, config, params)

temp_folder_path = config.root;
plots_folder_name = 'plots';
set(0,'DefaultFigureVisible','off');

if ~isfolder(fullfile(temp_folder_path, "snippets"))
    mkdir(fullfile(temp_folder_path, "snippets"));
else
    rmdir(fullfile(temp_folder_path, "snippets"),'s');
    mkdir(fullfile(temp_folder_path, "snippets"));
end

for fileID=1:size(segments.segments,2)
    segment = segments.segments(fileID).segment;
    name = segments.segments(fileID).name;
    signalStatus = segments.segments(fileID).signalStatus;
    % writing the metadata for file save at the end
    names{fileID,1} = name(1:end-4);

    [~, segment_m, indexes] = prepareInputForPlanner(segment);
    curveTypes = calculateCurveType(segment_m, indexes, name);
    frequencies = [linspace(0.016, 0.098, 5), linspace(0.114, 1, 10)];

    for freqID=1:length(frequencies)
        [P0_ma, ~, ~] = functional_calculateStraightLineOffset(segment_m, indexes, frequencies(freqID), mean(diff(segment_m(:,indexes.q_T0))));
        offsetError = -segment_m(:,indexes.c0)-P0_ma;
        thd = std(offsetError(curveTypes==1))*0.1;
        interactionPoints = abs(offsetError) > thd;
        interactionPoints(curveTypes ~=1) = false;
        interactionPoints = morphologyOpen(interactionPoints, 10);
        interactionPoints = morphologyClose(interactionPoints, 10);
        extrema = findExtremumPoints(offsetError, interactionPoints);

        %% cut snippets
        j = 1;
        k = 1;
        kpisNegOffs = [];
        kpisPosOffs = [];
        clear snippetsNegOffs snippetsPosOffs
    
        for i=1:size(extrema,1)
            snippetStart = extrema(i,1);
            % snippet start is modified to get the last point, where the thd
            % was exceeded before reaching the extremum
            offsetsBackwards = offsetError(snippetStart:-1:1);
            if (extrema(i,2) > 0)
                % this is a positive offset, searching for first point
                % backwards with negative sign
                thresholdExceed = find(offsetsBackwards <= 0, 1);
            else
                thresholdExceed = find(offsetsBackwards >= 0, 1);
            end
    
            snippetStop = snippetStart+find(abs(offsetError(snippetStart:end)) < 0.02,1); % the first point where c0 reaches its average after extremum is detected
            if (~isempty(thresholdExceed))
                snippetStart = snippetStart-thresholdExceed;
            end
            if (isempty(snippetStop))
               break;
            end
            compensationAborted = false;
            if (i<size(extrema,1))
                if (snippetStop > extrema(i+1,1))
                    % this is the situation, where something happened before
                    % getting back to the settle offset, therefore we simply
                    % neglect this compensation
                    compensationAborted = true;
                end
            end
            if (compensationAborted)
            else
                snippetStop = min(snippetStop, size(segment_m,1));
                snippetStart = max(snippetStart,1); % adding some init phase
                snippetStopExtended = min(snippetStop + (snippetStop-snippetStart)*2, size(segment_m,1));
                if (true)
                    if (extrema(i,2) > 0)
                        % the offset extremum is above the straight line offset
                        % of the driver
                        snippetsPosOffs(j).name = name;
                        snippetsPosOffs(j).snippet = segment_m(snippetStart:snippetStopExtended, :);
                        snippetsPosOffs(j).indexes = indexes;
                        snippetsPosOffs(j).offsetError = offsetError(snippetStart:snippetStopExtended);
                        snippetsPosOffs(j).P0_ma = P0_ma(snippetStart:snippetStopExtended);
                        snippetsPosOffs(j).relevantPoints.snippetStart = snippetStart;
                        snippetsPosOffs(j).relevantPoints.snippetStop = snippetStop;
                        snippetsPosOffs(j).relevantPoints.snippetStopExtended = snippetStopExtended;
                        snippetsPosOffs(j).relevantPoints.offsetExtremaLocation = extrema(i,1);
                        % for positive offset the following quantities are
                        % calculated&saved
                        % - minimum absolute relative lateral acceleration
                        % (right side acceleration)
                        % - minimum absolute relative lateral jerk
                        % (right side jerk)
                        % - minimum absolute curvature change (right side
                        % curvature gradient)
                        % - maximum lateral offset to the left side
                        kpisPosOffs(j, 1:5) = [min(segment_m(snippetStart:snippetStop, indexes.ayRel)) ...
                            min(segment_m(snippetStart:snippetStop, indexes.jyRel)) ...
                            min(segment_m(snippetStart:snippetStop, indexes.vehicleCurvatureChange)) ...
                            min(segment_m(snippetStart:snippetStop, indexes.TTCL)) ...
                            max(offsetError(snippetStart:snippetStop))];
                         j= j+1;
                    else
                        snippetsNegOffs(k).name = name;
                        snippetsNegOffs(k).snippet = segment_m(snippetStart:snippetStopExtended, :);
                        snippetsNegOffs(k).indexes = indexes;
                        snippetsNegOffs(k).offsetError = offsetError(snippetStart:snippetStopExtended);
                        snippetsNegOffs(k).P0_ma = P0_ma(snippetStart:snippetStopExtended);
                        snippetsNegOffs(k).relevantPoints.snippetStart = snippetStart;
                        snippetsNegOffs(k).relevantPoints.snippetStop = snippetStop;
                        snippetsNegOffs(k).relevantPoints.snippetStopExtended = snippetStopExtended;
                        snippetsNegOffs(k).relevantPoints.offsetExtremaLocation = extrema(i,1);
                        kpisNegOffs(k, 1:5) = [max(segment_m(snippetStart:snippetStop, indexes.ayRel)) ...
                            max(segment_m(snippetStart:snippetStop, indexes.jyRel)) ...
                            max(segment_m(snippetStart:snippetStop, indexes.vehicleCurvatureChange)) ...
                            min(segment_m(snippetStart:snippetStop, indexes.TTCL)) ...
                            min(offsetError(snippetStart:snippetStop))];
                        k = k+1;
                    end
                end
            end
        end 
        offsetErrorVector{freqID} = offsetError;
        pointsOfInteres{freqID} = zeros(size(offsetError,1),1);
        for i=1:length(snippetsNegOffs)
            if (i>1)
                compensationLength{freqID}(i) = (snippetsNegOffs(i).relevantPoints.snippetStop - max(snippetsNegOffs(i).relevantPoints.snippetStart, snippetsNegOffs(i-1).relevantPoints.snippetStop));
                pointsOfInteres{freqID}(max(snippetsNegOffs(i).relevantPoints.snippetStart, snippetsNegOffs(i-1).relevantPoints.snippetStop) : snippetsNegOffs(i).relevantPoints.snippetStop) = 1;
            else
                compensationLength{freqID}(i) = (snippetsNegOffs(i).relevantPoints.snippetStop - snippetsNegOffs(i).relevantPoints.snippetStart);
                pointsOfInteres{freqID}(snippetsNegOffs(i).relevantPoints.snippetStart:snippetsNegOffs(i).relevantPoints.snippetStop) = 1;
            end            
        end
        for i=length(snippetsNegOffs)+1:length(snippetsNegOffs)+length(snippetsPosOffs)
            if (i>length(snippetsNegOffs)+2)
                compensationLength{freqID}(i) = (snippetsPosOffs(i-length(snippetsNegOffs)).relevantPoints.snippetStop - max(snippetsPosOffs(i-length(snippetsNegOffs)-1).relevantPoints.snippetStop, snippetsPosOffs(i-length(snippetsNegOffs)).relevantPoints.snippetStart));
                pointsOfInteres{freqID}(max(snippetsPosOffs(i-length(snippetsNegOffs)-1).relevantPoints.snippetStop, snippetsPosOffs(i-length(snippetsNegOffs)).relevantPoints.snippetStart) : snippetsPosOffs(i-length(snippetsNegOffs)).relevantPoints.snippetStop) = 1;
            else
                compensationLength{freqID}(i) = (snippetsPosOffs(i-length(snippetsNegOffs)).relevantPoints.snippetStop - snippetsPosOffs(i-length(snippetsNegOffs)).relevantPoints.snippetStart);
                pointsOfInteres{freqID}(snippetsPosOffs(i-length(snippetsNegOffs)).relevantPoints.snippetStart : snippetsPosOffs(i-length(snippetsNegOffs)).relevantPoints.snippetStop) = 1;
            end
        end

        f = figure('units','normalized','outerposition',[0 0 1 1]);
        subplot(2,1,1);
        plot(segment_m(:,indexes.q_T0), -segment_m(:,indexes.c0), 'color', 'k');
        hold on;             grid on;
        plot(segment_m(:,indexes.q_T0), P0_ma, 'LineStyle','--', 'color', 'k');

        subplot(2,1,2);
        plot(segment_m(curveTypes==1,indexes.q_T0), offsetError(curveTypes==1), 'LineWidth',1, 'color', 'k');
        hold on; grid on;
        plot([0, segment_m(end, indexes.q_T0)], [thd, thd], 'color', 'g');
        plot([0, segment_m(end, indexes.q_T0)], [-thd, -thd], 'color', 'g');

        for i=1:length(snippetsNegOffs)
            
            plot(segment_m(snippetsNegOffs(i).relevantPoints.snippetStart:...
                snippetsNegOffs(i).relevantPoints.snippetStop,indexes.q_T0), ...
                offsetError(snippetsNegOffs(i).relevantPoints.snippetStart:...
                snippetsNegOffs(i).relevantPoints.snippetStop), ...
                'LineWidth', 1.5, 'LineStyle','--', 'color', 'b');

            plot(segment_m(snippetsNegOffs(i).relevantPoints.offsetExtremaLocation,indexes.q_T0), ...
                offsetError(snippetsNegOffs(i).relevantPoints.offsetExtremaLocation), ...
                'bo', 'LineWidth', 2);
        end

        for i=1:length(snippetsPosOffs)
            
            plot(segment_m(snippetsPosOffs(i).relevantPoints.snippetStart:...
                snippetsPosOffs(i).relevantPoints.snippetStop,indexes.q_T0), ...
                offsetError(snippetsPosOffs(i).relevantPoints.snippetStart:...
                snippetsPosOffs(i).relevantPoints.snippetStop), ...
                'LineWidth', 1.6, 'LineStyle',':', 'color', 'r');

            plot(segment_m(snippetsPosOffs(i).relevantPoints.offsetExtremaLocation,indexes.q_T0), ...
                offsetError(snippetsPosOffs(i).relevantPoints.offsetExtremaLocation), ...
                'ro', 'LineWidth', 2);
        end

        plot(segment_m(extrema(:,1),indexes.q_T0), ...
                offsetError(extrema(:,1)), ...
                'kx', 'LineWidth', 2);

        savefig(f, fullfile(temp_folder_path, plots_folder_name,...
                strcat('OffsetSnippets', snippetsNegOffs(1).name(1:end-4), '.fig')));
        saveas(f, fullfile(temp_folder_path, plots_folder_name,...
                        strcat('OffsetSnippets', snippetsNegOffs(1).name(1:end-4), '.png')));
        close(f);

    end

    f = figure();
    f.Position = [100 100 750 400];
    set(0,'defaultAxesFontSize',14);
    
    for i=1:length(compensationLength)
        subplot(1,3,1);
        [N, edges] = histcounts(compensationLength{i}, 'Normalization', 'probability');
        plot(0.5*(edges(1:end-1)+edges(2:end))*0.05, N, 'DisplayName',num2str(frequencies(i)));
        hold on; xlabel('time(s)'); grid on;
        title('histogram of snippet length');
        legend;
        data(i).lengths = [mean(compensationLength{i}), std(compensationLength{i})]*0.05;

        subplot(1,3,2);
        % ratio of length described by snippets
        plot(frequencies(i), numel(find(curveTypes==1 & pointsOfInteres{i}==1))/numel(find(curveTypes==1)), 'kx');
        grid on;
        hold on;
        title('ratio of length described by the snippets');
        data(i).coverage = numel(find(curveTypes==1 & pointsOfInteres{i}==1))/numel(find(curveTypes==1));

        subplot(1,3,3);
        [N, edges] = histcounts(offsetErrorVector{i}(pointsOfInteres{i}==1), 'Normalization', 'probability');
        plot(0.5*(edges(1:end-1)+edges(2:end)), N, 'DisplayName',num2str(frequencies(i)));
        hold on; xlabel('offset error(m)'); grid on;
        title('histogram of offset error');
        legend;
        data(i).offsetError = [mean(offsetErrorVector{i}(pointsOfInteres{i}==1)), std(offsetErrorVector{i}(pointsOfInteres{i}==1))];
    end
        
    savefig(f, fullfile(temp_folder_path, plots_folder_name,...
                    strcat('Window_Analysis', snippetsNegOffs(i).name(1:end-4),'_', num2str(i), '.fig')));
    saveas(f, fullfile(temp_folder_path, plots_folder_name,...
                    strcat('Window_Analysis', snippetsNegOffs(i).name(1:end-4),'_', num2str(i), '.png')));
    close(f);

    save(fullfile(temp_folder_path, plots_folder_name,...
                    strcat('Window_AnalysisData', snippetsNegOffs(i).name(1:end-4),'_', num2str(i), '.mat')), 'data');
end

end

function snippetModelByMpc(segments, config, params)
global snippetsNegOffs i interventionPoint snippetLocalTrajectory snippetLocalVehiclePath snippetLength parameters

temp_folder_path = config.root;
plots_folder_name = 'plots';
set(0,'DefaultFigureVisible','off');

if ~isfolder(fullfile(temp_folder_path, "snippets"))
    mkdir(fullfile(temp_folder_path, "snippets"));
else
    rmdir(fullfile(temp_folder_path, "snippets"),'s');
    mkdir(fullfile(temp_folder_path, "snippets"));
end

optParamFolder = dir(fullfile("C:\git\KDP_Igneczi\publikációk\LaneWandering\data\mpcOptimization\optimizedParametersLeftSuccessful\*.mat"));


for fileID=1:size(segments.segments,2)
    segment = segments.segments(fileID).segment;
    name = segments.segments(fileID).name;
    signalStatus = segments.segments(fileID).signalStatus;
    % writing the metadata for file save at the end
    names{fileID,1} = name(1:end-4);

    for i=1:length(optParamFolder)
        if (contains(optParamFolder(i).name, name(1:end-8)))
            optParams = load(fullfile(optParamFolder(i).folder, optParamFolder(i).name));
            optParams = optParams.OptimizedParameters;
            break;
        end
    end

    [~, segment_m, indexes] = prepareInputForPlanner(segment);
    curveTypes = calculateCurveType(segment_m, indexes, name);
    window = linspace(50, 500, 11); % with step size of 50 samples = 2.5s

    [P0_ma, P0, P0_fd, P0_oc, error_fd, error_oc, errorValues_fd, errorValues_oc] = calculateStraightLineOffset(segment_m, indexes, curveTypes, 180);
    %P0_ma = ones(size(segment_m,1),1)*mean(-segment_m(curveTypes==1, indexes.c0));
    %P0_ma = P0;
    offsetError = -segment_m(:,indexes.c0)-P0_ma;
    thd = std(offsetError(curveTypes==1))*0.5;
    interactionPoints = abs(offsetError) > thd;
    interactionPoints(curveTypes ~=1) = false;
    interactionPoints = morphologyOpen(interactionPoints, 10);
    interactionPoints = morphologyClose(interactionPoints, 10);
    extrema = findExtremumPoints(offsetError, interactionPoints);

    %% cut snippets
    j = 1;
    k = 1;
    kpisNegOffs = [];
    kpisPosOffs = [];
    %clear snippetsNegOffs snippetsPosOffs
    snippetsNegOffs = [];
    snippetsPosOffs = [];

    for i=1:size(extrema,1)
        snippetStart = extrema(i,1);
        % snippet start is modified to get the last point, where the thd
        % was exceeded before reaching the extremum
        offsetsBackwards = offsetError(snippetStart:-1:1);
        if (extrema(i,2) > 0)
            % this is a positive offset, searching for first point
            % backwards with negative sign
            thresholdExceed = find(offsetsBackwards <= 0, 1);
        else
            thresholdExceed = find(offsetsBackwards >= 0, 1);
        end

        snippetStop = snippetStart+find(abs(offsetError(snippetStart:end)) < 0.02,1); % the first point where c0 reaches its average after extremum is detected
        if (~isempty(thresholdExceed))
            snippetStart = snippetStart-thresholdExceed;
        end
        if (isempty(snippetStop))
           break;
        end
        compensationAborted = false;
        if (i<size(extrema,1))
            if (snippetStop > extrema(i+1,1))
                % this is the situation, where something happened before
                % getting back to the settle offset, therefore we simply
                % neglect this compensation
                compensationAborted = true;
            end
        end
        if (compensationAborted)
        else
            snippetStop = min(snippetStop, size(segment_m,1));
            snippetStart = max(snippetStart,1); % adding some init phase
            snippetStopExtended = min(snippetStop + (snippetStop-snippetStart)*2, size(segment_m,1));
            if ((snippetStop-snippetStart) <= 180)
                if (extrema(i,2) > 0)
                    % the offset extremum is above the straight line offset
                    % of the driver
                    snippetsPosOffs(j).name = name;
                    snippetsPosOffs(j).snippet = segment_m(snippetStart:snippetStopExtended, :);
                    snippetsPosOffs(j).indexes = indexes;
                    snippetsPosOffs(j).offsetError = offsetError(snippetStart:snippetStopExtended);
                    snippetsPosOffs(j).P0_ma = P0_ma(snippetStart:snippetStopExtended);
                    snippetsPosOffs(j).relevantPoints.snippetStart = snippetStart;
                    snippetsPosOffs(j).relevantPoints.snippetStop = snippetStop;
                    snippetsPosOffs(j).relevantPoints.snippetStopExtended = snippetStopExtended;
                    snippetsPosOffs(j).relevantPoints.offsetExtremaLocation = extrema(i,1);
                    % for positive offset the following quantities are
                    % calculated&saved
                    % - minimum absolute relative lateral acceleration
                    % (right side acceleration)
                    % - minimum absolute relative lateral jerk
                    % (right side jerk)
                    % - minimum absolute curvature change (right side
                    % curvature gradient)
                    % - maximum lateral offset to the left side
                    kpisPosOffs(j, 1:5) = [min(segment_m(snippetStart:snippetStop, indexes.ayRel)) ...
                        min(segment_m(snippetStart:snippetStop, indexes.jyRel)) ...
                        min(segment_m(snippetStart:snippetStop, indexes.vehicleCurvatureChange)) ...
                        min(segment_m(snippetStart:snippetStop, indexes.TTCL)) ...
                        max(offsetError(snippetStart:snippetStop))];
                     j= j+1;
                else
                    snippetsNegOffs(k).name = name;
                    snippetsNegOffs(k).snippet = segment_m(snippetStart:snippetStopExtended, :);
                    snippetsNegOffs(k).indexes = indexes;
                    snippetsNegOffs(k).offsetError = offsetError(snippetStart:snippetStopExtended);
                    snippetsNegOffs(k).P0_ma = P0_ma(snippetStart:snippetStopExtended);
                    snippetsNegOffs(k).relevantPoints.snippetStart = snippetStart;
                    snippetsNegOffs(k).relevantPoints.snippetStop = snippetStop;
                    snippetsNegOffs(k).relevantPoints.snippetStopExtended = snippetStopExtended;
                    snippetsNegOffs(k).relevantPoints.offsetExtremaLocation = extrema(i,1);
                    kpisNegOffs(k, 1:5) = [max(segment_m(snippetStart:snippetStop, indexes.ayRel)) ...
                        max(segment_m(snippetStart:snippetStop, indexes.jyRel)) ...
                        max(segment_m(snippetStart:snippetStop, indexes.vehicleCurvatureChange)) ...
                        min(segment_m(snippetStart:snippetStop, indexes.TTCL)) ...
                        min(offsetError(snippetStart:snippetStop))];
                    k = k+1;
                end
            end
        end
    end 

    % now rerunning the snippets with proper vehicle init state and calling mpc
    % afterwards to generate the wandering path for both the stray-away and the
    % compensate
    % prediction horizon is estimated based on compensation length
    compensationLength = 0;
    snippetsNegOffs = snippetsPosOffs;
    kpisNegOffs = kpisPosOffs;

    for i=1:length(snippetsNegOffs)
        compensationLength = compensationLength+(snippetsNegOffs(i).relevantPoints.snippetStop - snippetsNegOffs(i).relevantPoints.snippetStart);
    end

    %% resimulation
    parameters.Np = 45;
    parameters.Nc = 45;

    parameters.mode = 2; % 1 snippet - 1 calibration
    weightCombinations = [0 0 100; 0.01 0 100; 0.001 0 100];

    % mode=2: 1 snippet, 2 calibration phases
    weightCombinationsDrift = [0 0.01 100];
    weightCombinationsCompensation = [10 0 10];
    
    % mode=3: 1 snippet. 3 calibration phases
    weightCombinationsDriftBegin = [0 0.01 100; 0 0.01 10; 0 0.01 1];
    weightCombinationsIntervention = [10 0 10; 1 0 10; 0.1 0 10];
    weightCombinationsCompensationBegin = [0.01 0 10; 1 0 10; 0.1 0 10];

    % initial parameters for optimization
    x0 = [0; 100; 0; 100];

   for i=1:min(25,min(size(optParams,1), length(snippetsNegOffs)))

        dataSize = snippetsNegOffs(i).relevantPoints.snippetStop - snippetsNegOffs(i).relevantPoints.snippetStart;
        localOffsetIndex = max(1,snippetsNegOffs(i).relevantPoints.offsetExtremaLocation-snippetsNegOffs(i).relevantPoints.snippetStart);

        
        clear trajectory;
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

        for k=1:size(weightCombinationsDrift,1)
            % assinging dynamic parameters
            parameters.qDrift = weightCombinationsDrift(k,1:2);
            parameters.rDrift = weightCombinationsDrift(k,3);
            parameters.qCompensate = weightCombinationsCompensation(k,1:2);
            parameters.rCompensate = weightCombinationsCompensation(k,3);
            parameters.q = weightCombinations(k,1:2);
            parameters.rw = weightCombinations(k,3);
            parameters.qDriftBegin = weightCombinationsDriftBegin(k,1:2);
            parameters.rDriftBegin = weightCombinationsDriftBegin(k,3);
            parameters.qIntervention = weightCombinationsIntervention(k,1:2);
            parameters.rIntervention = weightCombinationsIntervention(k,3);
            parameters.qCompensationBegin = weightCombinationsCompensationBegin(k,1:2);
            parameters.rCompensationBegin = weightCombinationsCompensationBegin(k,3);

            options = optimset('TolX', 1.1e-6, 'TolFun', 1.1e-6); 
            options.OutputFcn = @custom_stop_fun;
            
            %[x, f] = fminsearch(@objectiveFcn, x0, [], [], [], [], zeros(4,1), [], [], options);
            % **************** OPTIMIZATION OF WEIGHTS ******************
            %[x, f] = fminsearch(@objectiveFcn, x0, options);

            %fprintf('Optimization of iteration %d has run\n',i);
            %fprintf('Final cost is %f, final parameters are [%f %f %f %f]\n', f, x);

            %x0 = x; % returning previous cycle value

%             parameters.qDrift = [x(1) 0];
%             parameters.rDrift = x(2);
%             parameters.qCompensate = [x(3) 0];
%             parameters.rCompensate = x(4);
% 
%             y = [x; f];
% 
%             save(fullfile(temp_folder_path, "snippets",...
%                             strcat('OptimizedParameters_', snippetsNegOffs(i).name(1:end-4),'_', num2str(i), '.mat')), "y");

            % ****************** RESIMULATION WITH SPECIFIC WEIGHTS*******
            % driver specific weight from external source
            parametersSpecific = table2array(readtable("C:\git\KDP_Igneczi\publikációk\LaneWandering\data\parameters.xlsx"));
%             for paramID=1:size(parametersSpecific,1)
%                 if (parametersSpecific(paramID,1) == str2num(name(3:5)))
%                     parameters.qDriftLeft = [parametersSpecific(paramID,2) 0];
%                     parameters.rDriftLeft = parametersSpecific(paramID,3);
%                     parameters.qCompensateLeft = [parametersSpecific(paramID,4) 0];
%                     parameters.rCompensateLeft = parametersSpecific(paramID,5);
%                     parameters.qDriftRight = [parametersSpecific(paramID,6) 0];
%                     parameters.rDriftRight = parametersSpecific(paramID,7);
%                     parameters.qCompensateRight = [parametersSpecific(paramID,8) 0];
%                     parameters.rCompensateRight = parametersSpecific(paramID,9);
%                     break;
%                 end
%             end     

            % **************** READ OF MOE/LOE parameters *****************
            load("C:\git\KDP_Igneczi\publikációk\LaneWandering\data\MOE.mat");
            load("C:\git\KDP_Igneczi\publikációk\LaneWandering\data\LOE.mat");

            
            parameters.qDrift = [optParams(i,1) 0];
            parameters.rDrift = optParams(i,2);
            parameters.qCompensate = [optParams(i,3) 0];
            parameters.rCompensate = optParams(i,4);
            
            [snippetsNegOffs(i).simulatedPos, interventionPointOE] = resimulateSnippet(snippetsNegOffs(i).snippet, snippetsNegOffs(i).indexes, ...
                interventionPoint, ...
                snippetLocalTrajectory, ...
                snippetLocalVehiclePath, ...
                snippetLength, ...
                parameters, ...
                MOE(fileID,:), ...
                LOE(fileID,:));

            % calculate the error between the planned and original paths
            ref(:,1) = snippetsNegOffs(i).simulatedPos(:,1);
            ref(:,2) = spline(snippetLocalVehiclePath(:,1), snippetLocalVehiclePath(:,2), ref(:,1));

            if (all(optParams(i,:)==0))
                nrms = 'nan';
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
        end

        % ******************* PLOTTING **********************************
        if (size(snippetsNegOffs(i).simulatedPos,1) > 5)

            f = figure();
            
            subplot(2,1,1);
            % path output
            plot(snippetsNegOffs(i).simulatedPos(:,1), snippetsNegOffs(i).simulatedPos(:,2), 'r', 'LineWidth', 2);
            hold on; grid on;
            plot(snippetLocalTrajectory(1:dataSize,1), snippetLocalTrajectory(1:dataSize,2), 'k');
            plot(snippetLocalVehiclePath(1:dataSize,1), snippetLocalVehiclePath(1:dataSize,2), 'b');
            plot([snippetLocalTrajectory(localOffsetIndex) snippetLocalTrajectory(localOffsetIndex)], ...
                [get(gca,'YLim')], 'k', 'LineWidth', 2);
            plot([snippetLocalVehiclePath(interventionPointOE,1) snippetLocalVehiclePath(interventionPointOE,1)], ...
                [get(gca,'YLim')], 'color', [192 192 192]/255, 'LineWidth', 2);
            

            lgd = legend('Simulation', 'Path', 'Human', 'Intervention Point', 'Intervention Point xOE');
            lgd.Orientation = 'vertical';
            lgd.Location = 'best';
            xlabel('x(m)', 'FontSize', 12); ylabel('y(m)', 'FontSize', 12);
            set(gca,'FontSize', 10);
            title('Path selection during drifting and compensation', 'FontSize',14);
    
            subplot(2,1,2);
            % offset output
            % simulated output resampled at X of reference
            pathY = spline(snippetsNegOffs(i).simulatedPos(:,1), snippetsNegOffs(i).simulatedPos(:,2), snippetLocalTrajectory(1:dataSize,1));
            simulatedOffset = pathY - snippetLocalTrajectory(1:dataSize,2);
            plot(snippetLocalTrajectory(1:dataSize,1), snippetsNegOffs(i).offsetError(1:dataSize), 'b');
            hold on; grid on;
            xlabel('x(m)', 'FontSize', 12);
            ylabel('offset error(m)', 'FontSize', 12);

            plot(snippetLocalTrajectory(1:dataSize,1), simulatedOffset, 'r');
            ylim([-0.5, 0.5]);
            plot([snippetLocalTrajectory(localOffsetIndex) snippetLocalTrajectory(localOffsetIndex)], ...
                [get(gca,'YLim')], 'k', 'LineWidth', 2);
            plot([snippetLocalTrajectory(interventionPointOE) snippetLocalTrajectory(interventionPointOE)], ...
                [get(gca,'YLim')], 'color', [192 192 192]/255, 'LineWidth', 2);
            
            set(gca,'FontSize', 10);
            lgd = legend('human', 'simulated');
            lgd.Orientation = 'horizontal';
            lgd.Location = 'south';
            title('Offset error', 'FontSize',14);

            savefig(f, fullfile(temp_folder_path, "snippets",...
                            strcat('Resimulated_snippets', snippetsNegOffs(i).name(1:end-4),'_', num2str(i), '.fig')));
            saveas(f, fullfile(temp_folder_path, "snippets",...
                            strcat('Resimulated_snippets', snippetsNegOffs(i).name(1:end-4),'_', num2str(i), '.png')));
            close(f);

            % *************** DIFFERENT PARAMETERS ************************

%             f = figure();
% 
%             for k=1:length(data)
%                 subplot(2,1,1);
%                 % path output
%                 plot(data(k).snippetsNegOffs.simulatedPos(:,1), data(k).snippetsNegOffs.simulatedPos(:,2), 'r', 'LineWidth', 2);
%                 hold on; grid on;
%                 plot(data(k).snippetLocalTrajectory(1:dataSize,1), data(k).snippetLocalTrajectory(1:dataSize,2), 'k');
%                 plot(data(k).snippetLocalVehiclePath(1:dataSize,1), data(k).snippetLocalVehiclePath(1:dataSize,2), 'b');
%                 plot([data(k).snippetLocalTrajectory(localOffsetIndex) data(k).snippetLocalTrajectory(localOffsetIndex)], ...
%                     [min(data(k).snippetLocalVehiclePath(1:dataSize,2)), max(data(k).snippetLocalVehiclePath(1:dataSize,2))], 'k', 'LineWidth', 2);
%                 lgd = legend('simulation', 'path', 'human', 'offsetExtremum');
%                 lgd.Orientation = 'horizontal';
%                 lgd.Location = 'south';
%                 xlabel('x(m)', 'FontSize', 12); ylabel('y(m)', 'FontSize', 12);
%                 set(gca,'FontSize', 10);
%                 title('Path selection during drifting and compensation', 'FontSize',14);
%             
%                 subplot(2,1,2);
%                 % offset output
%                 % simulated output resampled at X of reference
%                 pathY = spline(data(k).snippetsNegOffs.simulatedPos(:,1), data(k).snippetsNegOffs.simulatedPos(:,2), data(k).snippetLocalTrajectory(1:dataSize,1));
%                 simulatedOffset = pathY - data(k).snippetLocalTrajectory(1:dataSize,2);
%                 plot(data(k).snippetLocalTrajectory(1:dataSize,1), data(k).snippetsNegOffs.offsetError(1:dataSize), 'b');
%                 hold on; grid on;
%                 xlabel('x(m)', 'FontSize', 12);
%                 ylabel('offset error(m)', 'FontSize', 12);
%             
%                 plot(data(k).snippetLocalTrajectory(1:dataSize,1), simulatedOffset, 'r');
%                 plot([data(k).snippetLocalTrajectory(localOffsetIndex) data(k).snippetLocalTrajectory(localOffsetIndex)], ...
%                     [min(simulatedOffset), max(simulatedOffset)], 'k', 'LineWidth', 2);
%                 ylim([-0.25, 0.25]);
%                 set(gca,'FontSize', 10);
%                 lgd = legend('human', 'simulated');
%                 lgd.Orientation = 'horizontal';
%                 lgd.Location = 'south';
%                 title('Offset error', 'FontSize',14);
%             
%             end
%             
%             savefig(f, fullfile(temp_folder_path, "snippets",...
%                             strcat('Batch_Resimulated_snippets', snippetsNegOffs(i).name(1:end-4),'_', num2str(i), '.fig')));
%             saveas(f, fullfile(temp_folder_path, "snippets",...
%                             strcat('Batch_Resimulated_snippets', snippetsNegOffs(i).name(1:end-4),'_', num2str(i), '.png')));
%             close(f);
        end
   end
   save(fullfile(temp_folder_path, "snippets",...
                             strcat('Data', snippetsNegOffs(i).name(1:end-4), '.mat')), 'data');
   clear data;
end   

end

function stop = custom_stop_fun(~, optimValues, ~)
if optimValues.fval <=1e-4
    stop = true;
elseif optimValues.iteration > 80
        if optimValues.fval <=3e-4
            stop = true;
        else
            stop = false;
        end
elseif optimValues.iteration > 200
    stop = true;
else
    stop = false;
end
end

function f = objectiveFcn(x0)
global snippetsNegOffs i interventionPoint snippetLocalTrajectory snippetLocalVehiclePath snippetLength parameters

parameters.qDrift = [x0(1) 0];
            parameters.rDrift = x0(2);
            parameters.qCompensate = [x0(3) 0];
            parameters.rCompensate = x0(4);

    [simulatedPos] = resimulateSnippet(snippetsNegOffs(i).snippet, snippetsNegOffs(i).indexes, ...
                interventionPoint, ...
                snippetLocalTrajectory, ...
                snippetLocalVehiclePath, ...
                snippetLength, ...
                parameters);

    dataSize = size(snippetLocalVehiclePath,1);

    ref(:,1) = simulatedPos(:,1);
    ref(:,2) = spline(snippetLocalVehiclePath(:,1), snippetLocalVehiclePath(:,2), ref(:,1));

    f = sum((simulatedPos(:,2) - ref(:,2)).^2);
    f = f / dataSize;

    f = f + sum(-(min(x0,0))*100);

%     disp('current parameters:');
%     disp(x0);
%     disp('current loss:');
%     disp(f);
%     plot(simulatedPos(:,1), simulatedPos(:,2), ref(:,1), ref(:,2));
%     grid on;
%     title(num2str(f));
%     shg;
end

function [simulatedPos, interventionPoint] = resimulateSnippet(snippet, indexes, interventionPoint, snippetLocalTrajectory, snippetLocalVehiclePath, snippetEnd, parameters, MOE, LOE)
    
    localOffsetIndex = max(1,interventionPoint);

    vehicleState = initVehicleState(snippet, indexes, 1, snippetLocalVehiclePath);
    dT = mean(diff(snippet(:,indexes.q_T0)));
    
    x_k1 =  [0; 0; vehicleState.vx; 0; 0];
    j = 1;
    run = true;
    interventionPointFound = false;
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
        elseif (parameters.mode == 3)
            if (vehicleState.X < snippetLocalTrajectory(localOffsetIndex,1)*0.75)
                q = parameters.qDriftBegin;
                rw = parameters.rDriftBegin;
            elseif (vehicleState.X < snippetLocalTrajectory(localOffsetIndex,1)*1.25)
                q = parameters.qIntervention;
                rw = parameters.rIntervention;
            else
                q = parameters.qCompensationBegin;
                rw = parameters.rCompensationBegin;
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

function snippetStatisticsOffsetBased (segments, config, params)

for fileID=1:size(segments.segments,2)
    segment = segments.segments(fileID).segment;
    name = segments.segments(fileID).name;
    signalStatus = segments.segments(fileID).signalStatus;
    % writing the metadata for file save at the end
    names{fileID,1} = name(1:end-4);

    [~, segment_m, indexes] = prepareInputForPlanner(segment);
    curveTypes = calculateCurveType(segment_m, indexes, name);

    if (indexes.oncomingTraffic > 0)
        oncomingTraffic = movmean(segment_m(:,indexes.oncomingTraffic) > 0, 40);
    else
        oncomingTraffic = zeros(size(segment_m,1),1);
    end

    temp_folder_path = config.root;
    plots_folder_name = 'plots';
    set(0,'DefaultFigureVisible','off');

    [P0_ma, P0, P0_fd, P0_oc, error_fd, error_oc, errorValues_fd, errorValues_oc] = calculateStraightLineOffset(segment_m, indexes, curveTypes, 180);
    offsetError = -segment_m(:,indexes.c0)-P0_ma;
    thd = std(offsetError(curveTypes==1))*0.1;
    interactionPoints = abs(offsetError) > thd;
    interactionPoints(curveTypes ~=1) = false;
    interactionPoints = morphologyOpen(interactionPoints, 10);
    interactionPoints = morphologyClose(interactionPoints, 10);
    extrema = findExtremumPoints(offsetError, interactionPoints);

    indVar = -segment_m(:,indexes.c0);

    % cut snippets
    j = 1;
    k = 1;
    kpisNegOffs = [];
    kpisPosOffs = [];

    for i=1:size(extrema,1)
        snippetStart = extrema(i,1);
        % snippet start is modified to get the last point, where the thd
        % was exceeded before reaching the extremum
        offsetsBackwards = offsetError(snippetStart:-1:1);
        if (extrema(i,2) > 0)
            % this is a positive offset, searching for first point
            % backwards with negative sign
            thresholdExceed = find(offsetsBackwards <= 0, 1);
        else
            thresholdExceed = find(offsetsBackwards >= 0, 1);
        end

        snippetStop = snippetStart+find(abs(offsetError(snippetStart:end)) < 0.02,1); % the first point where c0 reaches its average after extremum is detected
        if (~isempty(thresholdExceed))
            snippetStart = snippetStart-thresholdExceed;
        end
        if (isempty(snippetStop))
           break;
        end
        compensationAborted = false;
        if (i<size(extrema,1))
            if (snippetStop > extrema(i+1,1))
                % this is the situation, where something happened before
                % getting back to the settle offset, therefore we simply
                % neglect this compensation
                compensationAborted = true;
            end
        end
        if (compensationAborted)
        else
            snippetStop = min(snippetStop, size(segment_m,1)); %adding some settling phase
            snippetStart = max(snippetStart,1); % adding some init phase
            if ((snippetStop-snippetStart) <= 180)
                if (extrema(i,2) > 0)
                    % the offset extremum is above the straight line offset
                    % of the driver
                    snippetsPosOffs(j).name = name;
                    snippetsPosOffs(j).snippet = segment_m(snippetStart:snippetStop, :);
                    snippetsPosOffs(j).indexes = indexes;
                    snippetsPosOffs(j).offsetError = offsetError(snippetStart:snippetStop);
                    % for positive offset the following quantities are
                    % calculated&saved
                    % - minimum absolute relative lateral acceleration
                    % (right side acceleration)
                    % - minimum absolute relative lateral jerk
                    % (right side jerk)
                    % - minimum absolute curvature change (right side
                    % curvature gradient)
                    % - maximum lateral offset to the left side
                    kpisPosOffs(j, 1:7) = [min(segment_m(snippetStart:snippetStop, indexes.ayRel)) ...
                        min(segment_m(snippetStart:extrema(i,1), indexes.jyRel)) ...
                        mean(segment_m(snippetStart:extrema(i,1), indexes.velocityX)) ...
                        max(offsetError(snippetStart:snippetStop)) ...
                        min(segment_m(snippetStart:extrema(i,1), indexes.veta)) ...
                        min(segment_m(extrema(i,1):snippetStop, indexes.jyRel)) ...
                        max(oncomingTraffic(snippetStart:snippetStop))];
                     j= j+1;
                else
                    snippetsNegOffs(k).name = name;
                    snippetsNegOffs(k).snippet = segment_m(snippetStart:snippetStop, :);
                    snippetsNegOffs(k).indexes = indexes;
                    snippetsNegOffs(k).offsetError = offsetError(snippetStart:snippetStop);
                    kpisNegOffs(k, 1:7) = [max(segment_m(snippetStart:snippetStop, indexes.ayRel)) ...
                        max(segment_m(snippetStart:extrema(i,1), indexes.jyRel)) ...
                        mean(segment_m(snippetStart:extrema(i,1), indexes.velocityX)) ...
                        min(offsetError(snippetStart:snippetStop)) ...
                        max(segment_m(snippetStart:extrema(i,1), indexes.veta)) ...
                        max(segment_m(extrema(i,1):snippetStop, indexes.jyRel)) ...
                        max(oncomingTraffic(snippetStart:snippetStop))];
                    k = k+1;
                end
            end
        end
    end   

    % TTCL vs. lane offset comparison plot
    f = figure();
    f.Position = [100, 100, 1050, 550];
    
    for i=1:length(kpisNegOffs)
        subplot(2,2,1);
        plot(snippetsNegOffs(i).snippet(:, snippetsNegOffs(i).indexes.q_T0) - snippetsNegOffs(i).snippet(1, snippetsNegOffs(i).indexes.q_T0), 1./snippetsNegOffs(i).snippet(:, snippetsNegOffs(i).indexes.TTCL));
        hold on; grid on;
        xlabel('time(s)'); ylabel('TTLC^-1(1/s)');
        subplot(2,2,3);
        plot(snippetsNegOffs(i).snippet(:, snippetsNegOffs(i).indexes.q_T0) - snippetsNegOffs(i).snippet(1, snippetsNegOffs(i).indexes.q_T0), snippetsNegOffs(i).offsetError);
        hold on; grid on;
        xlabel('time(s)'); ylabel('offset error(m)');
    end

    savefig(f, fullfile(temp_folder_path, plots_folder_name,...
                    strcat('TTCL_vs_offset_', name(1:end-4), '.fig')));
    saveas(f, fullfile(temp_folder_path, plots_folder_name,...
                    strcat('TTCL_vs_offset_', name(1:end-4), '.png')));
    close(f);

    % Correlation plots    
    f = figure();
    f.Position = [100, 100, 1050, 550];
    relevantIndex = 4;

    subplot(1,3,1);
    plot(kpisPosOffs(:,relevantIndex), kpisPosOffs(:,5), 'bo');
    xlabel(params.xlabelstr); ylabel("v_{y,Rel}(m/s^2)");
    grid on;
    hold on;
    plot(kpisNegOffs(:,relevantIndex), kpisNegOffs(:,5), 'bx');
    xlim([-0.8,0.8]); ylim([-0.5, 0.5]);
    title('v_{y,rel}');

    subplot(1,3,2);
    plot(kpisPosOffs(:,relevantIndex), kpisPosOffs(:,1), 'bo');
    xlabel(params.xlabelstr); ylabel("a_{y,Rel}(m/s^2)");
    grid on;
    hold on;
    plot(kpisNegOffs(:,relevantIndex), kpisNegOffs(:,1), 'bx');
    xlim([-0.8,0.8]); ylim([-0.5, 0.5]);
    title('a_{y,rel}');

    subplot(1,3,3);
    plot(kpisPosOffs(:,relevantIndex), kpisPosOffs(:,2), 'bo');
    xlabel(params.xlabelstr); ylabel("j_{y,Rel}(m/s^3)");
    title('j_{y,rel}');
    grid on; hold on;
    xlim([-0.5,0.5]); ylim([-0.3 0.3]);
    plot(kpisNegOffs(:,relevantIndex), kpisNegOffs(:,2), 'bx');

    savefig(f, fullfile(temp_folder_path, plots_folder_name,...
                    strcat('SnippetStatistics_', name(1:end-4), '.fig')));
    saveas(f, fullfile(temp_folder_path, plots_folder_name,...
                    strcat('SnippetStatistics_', name(1:end-4), '.png')));
    close(f);

    data(fileID).negOffsTTLC = kpisNegOffs(:,4);
    data(fileID).posOffsTTLC = kpisPosOffs(:,4);
    data(fileID).kpisNegOffs = kpisNegOffs;
    data(fileID).kpisPosOffs = kpisPosOffs;
    c = polyfit(kpisNegOffs(:,relevantIndex), kpisNegOffs(:,1),1);
    LOE(fileID,2) = c(1); 
    c = polyfit(kpisPosOffs(:,relevantIndex), kpisPosOffs(:,1),1);
    LOE(fileID,1) = c(1); 
end

save(fullfile(temp_folder_path, plots_folder_name,...
                    strcat('LOE.mat')), "LOE");

% data summary plots for all drivers
driversUnderTest = [1];
f = figure();
f.Position = [100, 100, 650, 650];
nPlotVertical = 3; nPlotHorizontal = 3;
for i=1:length(driversUnderTest)
    subplot(nPlotVertical, nPlotHorizontal, i);
    plot(data(driversUnderTest(i)).kpisNegOffs(:, relevantIndex), data(driversUnderTest(i)).kpisNegOffs(:, 5), 'bo');
    hold on;
    plot(data(driversUnderTest(i)).kpisPosOffs(:, relevantIndex), data(driversUnderTest(i)).kpisPosOffs(:, 5), 'bx');
    grid on;
    xlabel('offset error extrema(m)');
    ylabel('v_{y,rel}(m/s)');
    xlim([-0.5, 0.5]); ylim([-0.5, 0.5]);
    title(strcat('Driver-', num2str(driversUnderTest(i))));
end

savefig(f, fullfile(temp_folder_path, plots_folder_name,...
                    strcat('CollectedCorrelationVy.fig')));
saveas(f, fullfile(temp_folder_path, plots_folder_name,...
                strcat('CollectedCorrelationVy.png')));
close(f);

f = figure();
f.Position = [100, 100, 650, 650];
for i=1:length(driversUnderTest)
    subplot(nPlotVertical, nPlotHorizontal, i);
    plot(data(driversUnderTest(i)).kpisNegOffs(:, relevantIndex), data(driversUnderTest(i)).kpisNegOffs(:, 1), 'bo');
    hold on;
    plot(data(driversUnderTest(i)).kpisPosOffs(:, relevantIndex), data(driversUnderTest(i)).kpisPosOffs(:, 1), 'bx');
    grid on;
    xlabel('offset error extrema(m)');
    ylabel('a_{y,rel}(m/s^2)');
    xlim([-0.5, 0.5]); ylim([-0.5, 0.5]);
    title(strcat('Driver-', num2str(driversUnderTest(i))));
end

savefig(f, fullfile(temp_folder_path, plots_folder_name,...
                    strcat('CollectedCorrelationAy.fig')));
saveas(f, fullfile(temp_folder_path, plots_folder_name,...
                strcat('CollectedCorrelationAy.png')));
close(f);

f = figure();
f.Position = [100, 100, 650, 650];
for i=1:length(driversUnderTest)
    subplot(nPlotVertical, nPlotHorizontal, i);
    plot(data(driversUnderTest(i)).kpisNegOffs(:, relevantIndex), data(driversUnderTest(i)).kpisNegOffs(:, 2), 'bo');
    hold on;
    plot(data(driversUnderTest(i)).kpisPosOffs(:, relevantIndex), data(driversUnderTest(i)).kpisPosOffs(:, 2), 'bx');
    grid on;
    xlabel('offset error extrema(m)');
    ylabel('j_{y,rel}(m/s^2)');
    xlim([-0.5, 0.5]); ylim([-0.3, 0.3]);
    title(strcat('Driver-', num2str(driversUnderTest(i))));
end

savefig(f, fullfile(temp_folder_path, plots_folder_name,...
                    strcat('CollectedCorrelationJy.fig')));
saveas(f, fullfile(temp_folder_path, plots_folder_name,...
                strcat('CollectedCorrelationJy.png')));
close(f);

f = figure();
f.Position = [100, 100, 1080, 525];
for i=1:length(data)
    
    subplot(1,2,1);
    title('Distribution of left side compensation');
    [N, bins] = histcounts(data(i).posOffsTTLC);
    plot(0.5*(bins(1:end-1)+bins(2:end)), N/sum(N), 'DisplayName', strcat('Driver', num2str(i)));
    hold on; grid on;
    xlabel('offset error extrema(m)'); ylabel('Probability');
    xlim([0, 0.8]); ylim([0,1]);

    subplot(1,2,2);
    title('Distribution of right side compensation');
    [N, bins] = histcounts(data(i).negOffsTTLC);
    plot(0.5*(bins(1:end-1)+bins(2:end)), N/sum(N), 'DisplayName', strcat('Driver', num2str(i)));
    hold on; grid on;
    xlabel('offset error extrema(m)'); ylabel('Probability');
    xlim([-0.8, 0]); ylim([0,1]);

    % saving the mean c0 extremum
    MOE(i,1) = mean(data(i).posOffsTTLC);
    MOE(i,2) = mean(data(i).negOffsTTLC);
end

save(fullfile(temp_folder_path, plots_folder_name,...
                    strcat('MOE.mat')), "MOE");

savefig(f, fullfile(temp_folder_path, plots_folder_name,...
                    strcat('TTLC_Distribution.fig')));
    saveas(f, fullfile(temp_folder_path, plots_folder_name,...
                    strcat('TTLC_Distribution.png')));
    close(f);

end

function offsetPlotting(segments, config)
    for fileID=1:size(segments.segments,2)
        segment = segments.segments(fileID).segment;
        name = segments.segments(fileID).name;
        signalStatus = segments.segments(fileID).signalStatus;
        cutInfo = segments.segments(fileID).cutInfo;
        [~, segment_m, indexes] = prepareInputForPlanner(segment);
        [CLP, CLP_function, CLP_series] = calculateCLP(segment_m, indexes, 2, false);
    
        % straight line section evaluation
        curveTypes = calculateCurveType(segment_m, indexes, name);
    
        temp_folder_path = config.root;
        plots_folder_name = 'plots';
        set(0,'DefaultFigureVisible','off');
        f = figure();
        subplot(2,1,1);
        plot(-segment_m(curveTypes==1, indexes.c0));
        grid on;
        title('offset on straight lines'); ylim([-1,1]); ylabel('offset(m)');
        subplot(2,2,3);
        [N, edges] = histcounts(-segment_m(curveTypes==1, indexes.c0), 20);
        edges = (edges(1:end-1)+edges(2:end))/2;
        bar(edges,N/numel(find(curveTypes==1)));
        str = strcat('std:', num2str(std(-segment_m(curveTypes==1, indexes.c0))));
        text(0, 0.45,str, 'FontSize', 10);
        str = strcat('mean:', num2str(mean(-segment_m(curveTypes==1, indexes.c0))));
        text(0, 0.35,str, 'FontSize', 10);
        ylim([0,0.5]);
        xlim([-1,1]);

        subplot(2,2,4);
        % frequency characteristic of the curvature
        Y = fft(-segment_m(curveTypes==1, indexes.c0));
        P2 = abs(Y/numel(-segment_m(curveTypes==1, indexes.c0)));
        P1 = P2(1:numel(-segment_m(curveTypes==1, indexes.c0))/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f_ = 20*(0:(numel(-segment_m(curveTypes==1, indexes.c0))/2))/numel(-segment_m(curveTypes==1, indexes.c0));
        plot(f_,movmean(P1,5)); xlim([0,0.5]);
        title("Single-Sided Amplitude Spectrum of X(t)");
        xlabel("f (Hz)");
        ylabel("|P1(f)|");
    
        savefig(f, fullfile(temp_folder_path, plots_folder_name,...
                    strcat('StraightLineOffsets_', name(1:end-4), '.fig')));
        saveas(f, fullfile(temp_folder_path, plots_folder_name,...
                    strcat('StraightLineOffsets_', name(1:end-4), '.png')));
        close(f);

    end
end

function extrema = findExtremumPoints(y, inters)
    dy = diff(y);
    dy = [dy; dy(end)];
    posGrad = dy>0;
    negGrad = dy<0;
    maxima = diff(posGrad)<0;
    maxima = [maxima; false];
    minima = diff(negGrad)<0;
    minima = [minima; false];
    extrema = [minima|maxima y];
   % extrema are the global extrema. If inters is empty, this is returned.
   % if not, then inters are considered only.
   if (~isempty(inters))
       % 1 inter means a section of following points
       i = 1; j = 1;
       while i<=length(inters)
           if(inters(i)==1)
               interStop = find(inters(i:end)==0,1);
               interStop = i+interStop-1;
               if (~isempty(interStop))
                   if (mean(y(i:interStop-1)) < 0)
                       [ymin, xmin] = min(y(i:interStop-1));
                       extremaSimplified(j,1:2) = [i+xmin ymin];
                       j = j+1;
                   else
                       [ymax, xmax] = max(y(i:interStop-1));
                       extremaSimplified(j,1:2) = [i+xmax ymax];
                       j = j+1;
                   end
                   i = interStop;
               else
                   break;
               end
           else
               i = i+1;
           end
       end
       extrema = extremaSimplified;
   end
end

function curveTypes = calculateCurveType(segment_m, indexes, name)

global temp_folder_path plots_folder_name
    % return: curveTypes: 0 = unknown, 1=straight, 2=left, 3=right
    thd = 3.5e-4;
    curveTypes = zeros(size(segment_m,1),1);
    curvature = movmean(segment_m(:,indexes.c2)*2, 50);
    straightLine = abs(curvature) < thd;
    straightLine = morphologyOpen(straightLine, 100);
    straightLine = morphologyOpen(straightLine, 50);
    straightLine = morphologyClose(straightLine, 100);
    curveTypes(straightLine) = 1;
    
    leftCurve = curvature>=thd;
    leftCurve = morphologyClose(leftCurve,50);
    rightCurve = curvature<=-thd;
    rightCurve = morphologyClose(rightCurve,50);
    curveTransition = leftCurve&rightCurve;
    curveTypes(leftCurve&~curveTransition) = 2;
    curveTypes(rightCurve&~curveTransition) = 3;
    curveTypes(curveTransition) = 2.5;

    g = figure();

    plot(segment_m(:,indexes.X_abs), segment_m(:, indexes.Y_abs));
    grid on; hold on;
    plot(segment_m(curveTypes==1, indexes.X_abs), segment_m(curveTypes==1, indexes.Y_abs), 'ko');
    legend('original', 'straights');
    savefig(g, fullfile(temp_folder_path, plots_folder_name,...
                strcat('Map', name(1:end-4), '.fig')));
    saveas(g, fullfile(temp_folder_path, plots_folder_name,...
                strcat('Map_', name(1:end-4), '.png')));
    close(g);

end

function dataOut = morphologyOpen(dataIn, windowSize)
% this function is a morphology open filter with window size filter
dataOut = dataIn;
for i=1:length(dataIn)-windowSize
    if (~all(dataIn(i:i+windowSize)))
        dataOut(i:i+windowSize) = 0;
    end
end

end

function dataOut = morphologyClose(dataIn, windowSize)
% this function is a morphology close filter with window size filter
dataOut = dataIn;
for i=1:length(dataIn)-windowSize
    if (any(dataIn(i:i+windowSize)))
        dataOut(i:i+windowSize) = 1;
    end
end

end


