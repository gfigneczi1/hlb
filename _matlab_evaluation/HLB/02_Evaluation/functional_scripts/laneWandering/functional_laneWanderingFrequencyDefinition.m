function functional_laneWanderingFrequencyDefinition(segments, config)

temp_folder_path = config.root;
plots_folder_name = 'plots';
set(0,'DefaultFigureVisible','off');

if ~isfolder(fullfile(temp_folder_path, "plots", "frequencies"))
    mkdir(fullfile(temp_folder_path, "plots", "frequencies"));
else
    rmdir(fullfile(temp_folder_path, "plots", "frequencies"),'s');
    mkdir(fullfile(temp_folder_path, "plots", "frequencies"));
end

for fileID=1:size(segments.segments,2)
    segment = segments.segments(fileID).segment;
    name = segments.segments(fileID).name;
    signalStatus = segments.segments(fileID).signalStatus;
    % writing the metadata for file save at the end
    names{fileID,1} = name(1:end-4);

    [~, segment_m, indexes] = prepareInputForPlanner(segment);
    curveTypes = calculateCurveType(segment_m, indexes);
    frequencies = [linspace(0.016, 0.098, 5), linspace(0.114, 1, 10)];

    for freqID=1:length(frequencies)
        [P0_ma, offsetError, ~] = functional_calculateStraightLineOffset(segment_m, indexes, frequencies(freqID), mean(diff(segment_m(:,indexes.q_T0))));
        parameters.offsetErrorThd = std(offsetError(curveTypes==1))*0.5;
        parameters.id = name;
        parameters.indexes = indexes;

        % cut snippets
        [snippetsPosOffs, snippetsNegOffs, extrema] = functional_laneWanderingSnippetting(offsetError, P0_ma, segment_m, curveTypes==1, parameters);    
  
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
        
        set(f,'defaulttextInterpreter','latex') ;
        set(f, 'defaultAxesTickLabelInterpreter','latex');  
        set(f, 'defaultLegendInterpreter','latex');
        
        subplot(2,1,1);
        plot(segment_m(:,indexes.q_T0), -segment_m(:,indexes.c0), 'color', 'k');
        hold on;             grid on;
        plot(segment_m(:,indexes.q_T0), P0_ma, 'LineStyle','--', 'color', 'k');

        subplot(2,1,2);
        plot(segment_m(curveTypes==1,indexes.q_T0), offsetError(curveTypes==1), 'LineWidth',1, 'color', 'k');
        hold on; grid on;
        plot([0, segment_m(end, indexes.q_T0)], [parameters.offsetErrorThd, parameters.offsetErrorThd], 'color', 'g');
        plot([0, segment_m(end, indexes.q_T0)], [-parameters.offsetErrorThd, -parameters.offsetErrorThd], 'color', 'g');

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

        savefig(f, fullfile(temp_folder_path, plots_folder_name, "frequencies", ...
                strcat('OffsetSnippets', snippetsNegOffs(1).name(1:end-4), '.fig')));
        saveas(f, fullfile(temp_folder_path, plots_folder_name, "frequencies", ...
                        strcat('OffsetSnippets', snippetsNegOffs(1).name(1:end-4), '.png')));
        close(f);
        
        %% command line communication
        clc;
        fprintf("progress: file %d / %d\n", fileID,  size(segments.segments,2));
        fprintf("progress: frequency %d / %d\n", freqID,  length(frequencies));
        
    end

    %% Summary plot
    f = figure();
    f.Position = [100 100 400 950];
    set(0,'defaultAxesFontSize',14);
    set(f,'defaulttextInterpreter','latex') ;
    set(f, 'defaultAxesTickLabelInterpreter','latex');  
    set(f, 'defaultLegendInterpreter','latex');
    
    for i=1:length(compensationLength)
        subplot(3,1,1);
        [N, edges] = histcounts(compensationLength{i}, 'Normalization', 'probability');
        plot(0.5*(edges(1:end-1)+edges(2:end))*0.05, N, 'DisplayName',num2str(frequencies(i)));
        hold on; xlabel('Time(s)'); grid on;
        title('$\sigma_c$', 'FontSize', 14);
        legend('Orientation', 'horizontal', 'location', 'best');
        data(i).lengths = [mean(compensationLength{i}), std(compensationLength{i})]*0.05;

        subplot(3,1,2);
        % ratio of length described by snippets
        plot(frequencies(i), numel(find(curveTypes==1 & pointsOfInteres{i}==1))/numel(find(curveTypes==1)), 'kx');
        grid on;
        hold on;
        title('\textbf{Coverage}', 'FontSize', 14);
        data(i).coverage = numel(find(curveTypes==1 & pointsOfInteres{i}==1))/numel(find(curveTypes==1));

        subplot(3,1,3);
        [N, edges] = histcounts(offsetErrorVector{i}(pointsOfInteres{i}==1), 'Normalization', 'probability');
        plot(0.5*(edges(1:end-1)+edges(2:end)), N, 'DisplayName',num2str(frequencies(i)));
        hold on; xlabel('offset error(m)'); grid on;
        title('$\mu_{\Delta y_{error}}$', 'FontSize', 14);
        legend('Orientation', 'horizontal', 'location', 'best');
        data(i).offsetError = [mean(offsetErrorVector{i}(pointsOfInteres{i}==1)), std(offsetErrorVector{i}(pointsOfInteres{i}==1))];
    end
        
    savefig(f, fullfile(temp_folder_path, plots_folder_name, "frequencies", ...
                    strcat('Window_Analysis', snippetsNegOffs(i).name(1:end-4),'_', num2str(i), '.fig')));
    saveas(f, fullfile(temp_folder_path, plots_folder_name, "frequencies", ...
                    strcat('Window_Analysis', snippetsNegOffs(i).name(1:end-4),'_', num2str(i), '.png')));
    close(f);

    save(fullfile(temp_folder_path, plots_folder_name, "frequencies", ...
                    strcat('Window_AnalysisData', snippetsNegOffs(i).name(1:end-4),'_', num2str(i), '.mat')), 'data');
end

end

