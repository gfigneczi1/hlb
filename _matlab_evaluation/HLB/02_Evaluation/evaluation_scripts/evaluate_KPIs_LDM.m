function KPI = evaluate_KPIs_LDM(segment, driverModelOutput, name, config, curveID)

temp_folder_path = config.root;
plots_folder_name = 'plots';

% KPI 1: how close do we get to the lane edges / do we violate lane edges?
deviationFromLaneCenter = ((driverModelOutput.traj(:,1)-driverModelOutput.cor(:,1)).^2+(driverModelOutput.traj(:,2)-driverModelOutput.cor(:,2)).^2).^0.5;
deviationFromLaneCenterRef = ((driverModelOutput.ref(:,1)-driverModelOutput.cor(:,1)).^2+(driverModelOutput.ref(:,2)-driverModelOutput.cor(:,2)).^2).^0.5;
KPI.borderViolation = numel(find((deviationFromLaneCenter>(1.875-0.85)).*(abs(deviationFromLaneCenterRef) < 1.25))) / size(deviationFromLaneCenter,1)*100;
if (KPI.borderViolation == 0)
    KPI.borderDistance = (1.875-0.85) - max(deviationFromLaneCenter);
else
    KPI.borderDistance = 0;
end
% KPI 2: deviation compared to human trajectory
deviationFromHuman = ((driverModelOutput.traj(:,1)-driverModelOutput.ref(:,1)).^2+(driverModelOutput.traj(:,2)-driverModelOutput.ref(:,2)).^2).^0.5;
KPI.maxDeviationFromHuman = max(deviationFromHuman);
KPI.avgDeviationFromHuman = mean(deviationFromHuman);

% KPI 3: finding the curves only
% first cutting the entire data to subsection based on jumps in the
% corridor
curves = cutCurves(driverModelOutput);
KPI_report_array = [];
if (isempty(curveID))
    N = size(curves,2);
    N0 = 1;
else
    N = curveID;
    N0 = N;
end
for j = N0:N
        name_fig = strcat("TrajectoryDeviation_",num2str(j),"_",name,".fig");
        name_png = strcat("TrajectoryDeviation_",num2str(j),"_",name,".png");
        KPI.curves(j-N0+1).data = curves(j);
        KPI.curves(j-N0+1).curve = evaluation_driverModelQualitativeMeasuresCurves(curves(j),name_fig,name_png);
        fn = fieldnames(KPI.curves(j-N0+1).curve);
        KPI_report_new_row = [];
        for i = 1:length(fn)
            KPI_report_new_row = [KPI_report_new_row KPI.curves(j-N0+1).curve.(fn{i})];
        end
        KPI_report_array = [KPI_report_array; KPI_report_new_row];
end
curveTraj = []; curveCor = []; curveRef=[];
for i=1:length(KPI.curves)
    curveTraj = [curveTraj; KPI.curves(i).data.traj];
    curveCor = [curveCor; KPI.curves(i).data.cor];
    curveRef = [curveRef; KPI.curves(i).data.ref];
end

plot_driverModelTrajectoryHistogram([],curveTraj,curveRef, strcat(name,"_deviationHistogramHumanInCurve",num2str(curveID)));

%% Save the mat file and the csv report
KPI_report_csv = array2table(KPI_report_array,'VariableNames',fn);
if (isempty(curveID))
    writetable(KPI_report_csv,fullfile(temp_folder_path, plots_folder_name,strcat(name,'.xlsx')));
else
    writetable(KPI_report_csv,fullfile(temp_folder_path, plots_folder_name,strcat(name,num2str(curveID),'.xlsx')));
end
set(0,'DefaultFigureVisible','on');
save(fullfile(temp_folder_path, plots_folder_name,strcat(name,num2str(curveID),".mat")),'KPI');

end

