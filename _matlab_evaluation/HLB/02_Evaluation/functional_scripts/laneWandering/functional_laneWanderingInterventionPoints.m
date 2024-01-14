function functional_laneWanderingInterventionPoints(segments, config)

for fileID=1:size(segments.segments,2)
    segment = segments.segments(fileID).segment;
    name = segments.segments(fileID).name;
    signalStatus = segments.segments(fileID).signalStatus;
    % writing the metadata for file save at the end
    names{fileID,1} = name(1:end-4);

    [~, segment_m, indexes] = prepareInputForPlanner(segment);
    curveTypes = calculateCurveType(segment_m, indexes);

    temp_folder_path = config.root;
    plots_folder_name = 'plots';
    set(0,'DefaultFigureVisible','off');

    [P0_ma, offsetError, ~] = functional_calculateStraightLineOffset(segment_m, indexes, 0.111, mean(diff(segment_m(:,indexes.q_T0))));
    
    parameters.offsetErrorThd = std(offsetError(curveTypes==1))*0.5;
    parameters.id = name;
    parameters.indexes = indexes;
    
    [~, ~, ~, kpisPosOffs, kpisNegOffs] = functional_laneWanderingSnippetting(offsetError, P0_ma, segment_m, curveTypes==1, parameters);   

    data(fileID).negOffsTTLC = kpisNegOffs(:,5); % 4 - TTCL, 5 - offset error
    data(fileID).posOffsTTLC = kpisPosOffs(:,5);
    data(fileID).kpisNegOffs = kpisNegOffs;
    data(fileID).kpisPosOffs = kpisPosOffs;
    
    %% command line communication
    clc;
    fprintf('Intervention Point Calculation\n');
    fprintf('Files %d / %d', fileID, size(segments.segments,2));    
end

f = figure();
f.Position = [100, 100, 525, 700];

set(f,'defaulttextInterpreter','latex') ;
set(f, 'defaultAxesTickLabelInterpreter','latex');  
set(f, 'defaultLegendInterpreter','latex');

for i=1:length(data)    
    subplot(2,1,1);
    [N, bins] = histcounts(data(i).posOffsTTLC);
    plot(0.5*(bins(1:end-1)+bins(2:end)), N/sum(N), 'DisplayName', strcat('Driver', num2str(i)));
    hold on; grid on;
    xlabel('$A^{dr}(m)$'); ylabel('$p(A^{dr})$');
    xlim([0, 0.8]); ylim([0,1]);
    legend('Location', 'best','NumColumns',2);
    set(gca, 'FontSize', 12);
    title('\textbf{Distribution of left side compensation}', 'FontSize', 14);
    xticks([0:0.1:0.8]);

    subplot(2,1,2);    
    [N, bins] = histcounts(data(i).negOffsTTLC);
    plot(0.5*(bins(1:end-1)+bins(2:end)), N/sum(N), 'DisplayName', strcat('Driver', num2str(i)));
    hold on; grid on;
    xlabel('$A^{dr}(m)$'); ylabel('$p(A^{dr})$');
    xlim([-0.8, 0]); ylim([0,1]);
    legend('Location', 'best','NumColumns',2);
    set(gca, 'FontSize', 12);
    title('\textbf{Distribution of right side compensation}', 'FontSize', 14);
    xticks([-0.8:0.1:0]);
end

savefig(f, fullfile(temp_folder_path, plots_folder_name,...
                    strcat('InterventionPoint_Distribution.fig')));
    saveas(f, fullfile(temp_folder_path, plots_folder_name,...
                    strcat('InterventionPoint_Distribution.png')));
    close(f);

end

