function plot_nodePointDefinition(outputStructs, options, config)
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
for i=1:numberOfEvals
    outputStruct = outputStructs{i};
    % create the folder for the eval
    plot_folder_path = fullfile(root_folder_path, strcat('Eval_',num2str(i)));
    mkdir(plot_folder_path);
    
    % 1. Surfaces of the parameter dependency
    if (options.parameterDependency == "true")
        f = figure('DefaultAxesFontSize',14);
        surf(outputStruct.paramDependency.P2, outputStruct.paramDependency.P3, outputStruct.paramDependency.surf_1);
        xlabel('P2(m)'); ylabel('P3(m)'); zlabel('Cost value');
        saveas(f, fullfile(plot_folder_path, 'P2_P3_dependency.png'));
        savefig(f,fullfile(plot_folder_path, 'P2_P3_dependency.fig'));
        close(f);
        f = figure('DefaultAxesFontSize',14);
        surf(outputStruct.paramDependency.P1, outputStruct.paramDependency.P3, outputStruct.paramDependency.surf_2);
        xlabel('P1(m)'); ylabel('P3(m)'); zlabel('Cost value');
        saveas(f, fullfile(plot_folder_path, 'P1_P3_dependency.png'));
        savefig(f,fullfile(plot_folder_path, 'P1_P3_dependency.fig'));
        close(f);
        f = figure('DefaultAxesFontSize',14);
        surf(outputStruct.paramDependency.P1, outputStruct.paramDependency.P2, outputStruct.paramDependency.surf_3);
        xlabel('P1(m)'); ylabel('P2(m)'); zlabel('Cost value');
        saveas(f, fullfile(plot_folder_path, 'P1_P2_dependency.png'));
        savefig(f,fullfile(plot_folder_path, 'P1_P2_dependency.fig'));
        close(f);
    end
    
    % 2. Init value dependency of the summed cost
    f = figure('DefaultAxesFontSize',14);
    plot(outputStruct.initvalueCosts(:,1),'-s','MarkerSize',10,'LineWidth',2);
    xlim([1, options.independentInitialValues]);
    grid on;
    xlabel('Init parameter index'); ylabel('Cost'); title('Cost value with multiple init values');
    saveas(f, fullfile(plot_folder_path, 'init_dependency.png'));
    savefig(f, fullfile(plot_folder_path, 'init_dependency.fig'));
    close(f);
    
    % 3. Driver's results in terms of final optimized parameter vectors and
    % their accuracy (mean, std and maximum) - no plot, only a mat file
    deviation = ((outputStruct.traj(:,2) - outputStruct.ref(:,2)).^2 ...
        + (outputStruct.traj(:,1) - outputStruct.ref(:,1)).^2).^0.5;
    driver_report.mean = mean(deviation);
    driver_report.std = std(deviation);
    driver_report.max = max(deviation);
    driver_report.f = outputStruct.f;
    driver_report.P3 = max(10,outputStruct.P3in*250);
    x3(1)= driver_report.P3(1);
    for j=2:length(driver_report.P3)
        x3 = [x3 x3(j-1)+driver_report.P3(j)];
    end
    driver_report.x3 = x3;
    clear x3;
    save(fullfile(plot_folder_path,'driver_report.mat'),'driver_report');
    
    disp(strcat('The statistical results are:'));
    disp(strcat('mean=', num2str(driver_report.mean)));
    disp(strcat('std=', num2str(driver_report.std)));
    disp(strcat('max=', num2str(driver_report.max)));

end



end

