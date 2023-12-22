function plot_driverModelLearning(driverModelOutput,KPI, config, segment, name, curveID)

temp_folder_path = config.root;
plots_folder_name = 'plots';

% Plotting histogram of lane center deviation
plot_driverModelTrajectoryHistogram(segment,driverModelOutput.traj,driverModelOutput.cor, strcat(name,"_deviationHistogramLaneCenter"));

% Plotting histogram of human trajectory deviation
plot_driverModelTrajectoryHistogram(segment,driverModelOutput.traj,driverModelOutput.ref, strcat(name,"_deviationHistogramHuman")); 

%% CURVE related plots based on KPI.curves
set(0,'DefaultFigureVisible','off');

if (isempty(curveID))
    N = size(KPI.curves,2);
    N0 = 1;
else
    N = curveID;
    N0 = N;
end

for i=N0:N
    name_fig = strcat("TrajectoryDeviation_",num2str(i),"_",name,".fig");
    name_png = strcat("TrajectoryDeviation_",num2str(i),"_",name,".png");
    
    curve = KPI.curves(i-N0+1).data;
    KPI_results = KPI.curves(i-N0+1).curve;
    % Data calculation
    T = [cos(curve.orient(1)) sin(curve.orient(1)); -sin(curve.orient(1)) cos(curve.orient(1))];
    ref = [curve.ref(:,1)-curve.cor(1,1) curve.ref(:,2)-curve.cor(1,2)]*T';
    cor = [curve.cor(:,1)-curve.cor(1,1) curve.cor(:,2)-curve.cor(1,2)]*T';
    traj = [curve.traj(:,1)-curve.cor(1,1) curve.traj(:,2)-curve.cor(1,2)]*T';
    orient = curve.orient - curve.orient(1);

    refSignedDeviation = sign(ref(:,2)-cor(:,2)).*(((ref(:,1)-cor(:,1)).^2+(ref(:,2)-cor(:,2)).^2).^0.5);
    trajSignedDeviation = sign(traj(:,2)-cor(:,2)).*(((traj(:,1)-cor(:,1)).^2+(traj(:,2)-cor(:,2)).^2).^0.5);

    refCurvature = -1*movmean([0; movmean([0; diff(diff(curve.ref(:,2))./diff(curve.ref(:,1)))]./diff(curve.ref(:,1)),250)], 50);
    trajCurvature = -1*movmean([0; movmean([0; diff(diff(curve.traj(:,2))./diff(curve.traj(:,1)))]./diff(curve.traj(:,1)),250)], 50);
    corCurvature = -1*movmean([0; movmean([0; diff(diff(curve.cor(:,2))./diff(curve.cor(:,1)))]./diff(curve.cor(:,1)),250)],50);

    corLocal = (curve.cor-curve.cor(1,:))*T';
    corLeftLocal = (curve.corLeft-curve.cor(1,:))*T';
    corRightLocal = (curve.corRight-curve.cor(1,:))*T';
    
    % Adding deviation plot (together with trajectory) - MAIN Plot
    deviationPlot = figure('Name', 'Trajectory deviation', 'NumberTitle','off');
    subplot(2,1,1);
    p = plot(corLocal(:,1), refSignedDeviation, corLocal(:,1), trajSignedDeviation, [corLocal(1,1) corLocal(end,1)],[0 0]);
    set(p(1),'LineWidth',1.5,'color','b');
    set(p(2),'LineWidth',1.5,'color','r');
    set(p(3),'LineWidth',1,'color',[17 17 17]/255);
    grid on; title('Lateral offset to midlane','FontSize',14);
    xlabel('X-coordinates(m)','FontSize',12); ylabel('Lateral offset(m)','FontSize',12);
    legend('ref to cor','traj to cor','FontSize',10);
    xlim([0 corLocal(end,1)]);
    ylim([-1.25 1.25]);
    xticks([0:200:corLocal(end,1)]); yticks([-1.25:0.25:1.25]);
    subplot(2,1,2);

    p = plot(corLocal(:,1),corLocal(:,2),corLeftLocal(:,1),corLeftLocal(:,2),corRightLocal(:,1),corRightLocal(:,2));
    grid on; axis equal;
    set(p(1),'LineWidth',1.5,'color','k');
    set(p(2),'LineWidth',1.5,'color',[17 17 17]/255,'LineStyle','--');
    set(p(3),'LineWidth',1.5,'color',[17 17 17]/255,'LineStyle','--');
    xlabel('X-coordinates (m)','FontSize',12);
    ylabel('Y-coordinates (m)','FontSize',12);
    title('Curve trajectory','FontSize',14);
    xticks([0:200:corLocal(end,1)]);

    annotation('textbox', [0.15, 0.85, 0.4, 0.05],...
                'String', strcat("Side correctness: ",num2str(KPI_results.sideCorrectness*100),"%"), ...
                'FontSize',10);

    savefig(deviationPlot, fullfile(temp_folder_path, plots_folder_name,...
            name_fig));
    saveas(deviationPlot, fullfile(temp_folder_path, plots_folder_name,...
           name_png));

    deviationPlot = figure('Name', 'Curvature deviation', 'NumberTitle','off');
    subplot(2,1,1);
    p = plot(corLocal(:,1), refCurvature*1000 , corLocal(:,1), trajCurvature*1000, corLocal(:,1), corCurvature*1000, [corLocal(1,1) corLocal(end,1)],[0 0]);
    set(p(1),'LineWidth',1.5,'color','b');
    set(p(2),'LineWidth',1.5,'color','r');
    set(p(3),'LineWidth',1.5,'color','k','LineStyle','--');
    set(p(4),'LineWidth',1,'color',[17 17 17]/255);
    xlim([0 corLocal(end,1)]);
    ylim([-4 4]);
    grid on;
    xlabel('X-coordinates (m)','FontSize',12); ylabel('Normalized curvature (1/m)','FontSize',12);
    title("Curvature of different trajectories");
    legend('ref','traj', 'cor');

    subplot(2,1,2);
    p = plot(corLocal(:,1), (refCurvature - corCurvature)*1000, corLocal(:,1), (trajCurvature-corCurvature)*1000);
    set(p(1),'LineWidth',1.5,'color','b');
    set(p(2),'LineWidth',1.5,'color','r');
    grid on;
    xlim([0 corLocal(end,1)]);
    ylim([-1 1]);
    xlabel('X-coordinates (m)','FontSize',12); ylabel('Normalized curvature difference (1/m)','FontSize',12);
    title("Curvature difference to corridor");
    legend('ref to cor','traj to cor');

    name_fig = strcat('curvature',name_fig);
    name_png = strcat('curvature',name_png);

    savefig(deviationPlot, fullfile(temp_folder_path, plots_folder_name,...
            name_fig));
    saveas(deviationPlot, fullfile(temp_folder_path, plots_folder_name,...
            name_png));
end
set(0,'DefaultFigureVisible','on');


end

