function plot_nodePointNumberAnalysis(outputStructs, options, config)
% The following plots shall be created:
% 1. Surfaces of the parameter dependency of the cost value (node point
% distance parameter vector)
% 2. Init value dependency of the summed cost
% 3. Driver's results in terms of final optimized parameter vectors and
% their accuracy (mean, std and maximum) - this driver report
% - we save metadata regarding evaluation
% 4. Driver's results in terms of common parameter vector and
% their accuracy (mean, std and maximum)

% INPUTS of the plots:
% 1. Planned global trajectory X-Y points (traj)
% 2. Global reference trajectory which is the driven path of the human
% driver (ref)
% 3. f_sum_array: an NxN array and the parameter vectors attached to it,
% which are 1xN sized 
% 4. P vector used during the curve fitting (used in the cost function)

% @copyright (C) 2022 Robert Bosch GmbH.
% The reproduction, distribution and utilization of this file as
% well as the communication of its contents to others without express
% authorization is prohibited. Offenders will be held liable for the
% payment of damages. All rights reserved in the event of the grant
% of a patent, utility model or design.
% @version 1.0

set(0,'DefaultFigureVisible','off');
root_folder_path = fullfile(config.root,'plots');
numberOfEvals = size(outputStructs,2);
markers = {'o', '+', '*', 'x', '.'};

inaccuracy = 0.03; % in meters, from localization and lane edge detection

f = figure();
f.Position = [100, 100, 640, 320];

set(f,'defaulttextInterpreter','latex') ;
set(f, 'defaultAxesTickLabelInterpreter','latex');  
set(f, 'defaultLegendInterpreter','latex');

if(~exist("root_folder_path"))
    mkdir(root_folder_path);
end

for i=1:numberOfEvals
    outputStruct = outputStructs{i};

    
    % plot the median and mean NRMS values of drivers, for various node
    % points
    meanNRMS = mean(outputStruct.nodePointCheckFvals');
    medianNRMS = median(outputStruct.nodePointCheckFvals');
   
    plot(meanNRMS, 'Marker',markers{i}, 'color', 'k', 'LineWidth', 1, 'DisplayName', append("Driver", num2str(i))); hold on; grid on;

    xlabel('$n_{nodePoints}$'); ylabel('$RMS_{\Delta y}(m)$'); 
    set(gca,'FontSize', 14);
    ylim([0,0.25]);
    xlim([1, options.parameters.Np]);
    xticks(1:1:size(outputStruct.nodePointCheckFvals,2));
    title('\textbf{$RMS_{\Delta y}$ for multiple node point selection}', 'FontSize', 15.4);
end

h = fill([1, options.parameters.Np options.parameters.Np 1], [0, 0, inaccuracy, inaccuracy], 'yellow', 'FaceAlpha',0.3, 'EdgeColor','none', 'DisplayName', 'Offset 95% inaccuracy zone');

legend('Location', 'northeast');

saveas(f, fullfile(root_folder_path, 'nodePointAnalysis.png'));
savefig(f,fullfile(root_folder_path, 'nodePointAnalysis.fig'));
close(f);
end

