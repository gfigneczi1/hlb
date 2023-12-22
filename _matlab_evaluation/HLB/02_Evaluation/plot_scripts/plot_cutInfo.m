function plot_cutInfo(config, rawData, cutInfo, token)
    temp_folder_path = config.root;
    plots_folder_name = 'plots';
    set(0,'DefaultFigureVisible','off');
    
    X_abs = rawData.LongPos_abs * 40075000 .* cos(rawData.LatPos_abs*pi()/180) / 360;
    Y_abs = rawData.LatPos_abs * 111.32*1000;
    
    f = figure();
    
    p = plot(X_abs, Y_abs);
    hold on;
    fn = fieldnames(cutInfo);
    all_marks = {'bo', 'bx', 'k*', 'g+'};
    l{1} = 'route';
    j = 2;
    for i=1:length(fn)
        if (numel(find(cutInfo.(fn{i})==1)) > 0)
            plot(X_abs(cutInfo.(fn{i})==1),Y_abs(cutInfo.(fn{i})==1), all_marks{i});
            l{j} = fn{i};
            j=j+1;
        end
    end
    axis equal;
    grid on; 
    legend(l, 'location','best');
    
    savefig(f, fullfile(temp_folder_path, plots_folder_name,...
        strcat('TrajectoryComparisonCutInfo',token,'.fig')));
    saveas(f, fullfile(temp_folder_path, plots_folder_name,...
        strcat('TrajectoryComparisonCutInfo',token,'.png')));
    
    close(f);
    set(0,'DefaultFigureVisible','on');
end

