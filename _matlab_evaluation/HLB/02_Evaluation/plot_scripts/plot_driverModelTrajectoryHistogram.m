function plot_driverModelTrajectoryHistogram(segment, a, b, name)
%DRIVERMODELTRAJECTORYHISTOGRAM Summary of this function goes here
%   Detailed explanation goes here

temp_folder_path = fullfile('..','..','_temp');
plots_folder_name = 'plots';
set(0,'DefaultFigureVisible','off');

meas_X = b(:,1);
meas_Y = b(:,2);
planned_X = a(:,1);
planned_Y = a(:,2);
orientation = diff(meas_Y)./diff(meas_X);
orientation = [orientation; orientation(end)];
orientation = atan(orientation);
indeces = find(isnan(orientation));
if (~isempty(indeces))
    orientation(find(isnan(orientation))) = 0;
    diffs = (meas_Y-planned_Y)./cos(orientation);
    diffs(indeces) = 0;
else
    diffs = (meas_Y-planned_Y)./cos(orientation);
end
%diffs = ((meas_X-planned_X).^2+(meas_Y-planned_Y).^2).^0.5;

std_calc = std(diffs);
rms_calc = rms(diffs);
mean_calc = mean(diffs);

histPlot = figure('Name', 'Relative trajectory deviation histogram', 'NumberTitle','off');
h = histogram(diffs, 'DisplayName', 'Deviation', 'Normalization', 'Probability', 'FaceAlpha',1);
patch([mean_calc-std_calc mean_calc+std_calc...
       mean_calc+std_calc mean_calc-std_calc], [max(ylim) max(ylim) 0 0],...
      'r', 'DisplayName',sprintf('STD: %.2f', std_calc), 'FaceAlpha', .5,...
      'edgecolor', 'none' );
title('Relative trajectory deviation histogram');
ylabel('Probability');
xlabel('Error [m]');
hold on; grid on;
xline(rms_calc,'--g', 'LineWidth', 2.5, 'DisplayName', sprintf('RMS: %.2f', rms_calc));
xline(mean_calc,'--y', 'LineWidth', 2.5, 'DisplayName', sprintf('Mean: %.2f', mean_calc));
uistack(h, 'top');
legend('show');
hold off;
name_fig = strcat(name,'.fig');
name_png = strcat(name,'.png');
set(histPlot, 'CreateFcn', 'set(gcbo,''Visible'',''on'')'); 
savefig(histPlot, fullfile(temp_folder_path, plots_folder_name,...
        name_fig));
saveas(histPlot, fullfile(temp_folder_path, plots_folder_name,...
       name_png));
set(0,'DefaultFigureVisible','on');
end

