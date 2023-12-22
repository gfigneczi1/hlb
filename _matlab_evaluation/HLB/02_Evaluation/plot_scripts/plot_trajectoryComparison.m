function plot_trajectoryComparison(config,ref,traj,cor,corLeft,corRight, token)
    temp_folder_path = config.root;
    plots_folder_name = 'plots';
    set(0,'DefaultFigureVisible','off');
    
    f = figure();
    
    p = plot(ref(:,1),ref(:,2),cor(:,1),cor(:,2),corLeft(:,1),corLeft(:,2),corRight(:,1),corRight(:,2), traj(:,1),traj(:,2));
    grid on; 
    xlabel('Global coordinates X (m)');
    ylabel('Global coordinates Y (m)');
    legend('Reference','MidLane','Left','Right','Trajectory');
    set(p(1),'LineWidth',2,'color','b','LineStyle','--');
    set(p(2),'LineWidth',2.5,'color','k','LineStyle','--');
    set(p(3),'LineWidth',1,'color','g','LineStyle',':');
    set(p(4),'LineWidth',1,'color','g','LineStyle',':');
    set(p(5),'LineWidth',1,'color','c');
    
    savefig(f, fullfile(temp_folder_path, plots_folder_name,...
        strcat('TrajectoryComparison',token,'.fig')));
    saveas(f, fullfile(temp_folder_path, plots_folder_name,...
        strcat('TrajectoryComparison',token,'.png')));
    
    close(f);
    set(0,'DefaultFigureVisible','on');
end

