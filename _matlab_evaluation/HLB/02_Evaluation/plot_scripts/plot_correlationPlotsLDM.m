function plot_correlationPlotsLDM(config,GT_U, GT_Y, token)

temp_folder_path = config.root;
plots_folder_name = 'plots';
set(0,'DefaultFigureVisible','off');

disp(strcat('Average straight line offset for driver', token));
disp(num2str(mean(GT_Y(2,abs(GT_U(1,:)) < 0.25))));

for i=1:size(GT_U,1)
    for j=1:size(GT_Y,1)-1 %-1 offset as first node point is always the origo
        f = figure();
        plot(GT_U(i,(abs(GT_U(i,:))>0.5)),GT_Y(j+1,(abs(GT_U(i,:))>0.5)),'kx','MarkerSize',3); 
        hold on;
        if numel(find(abs(GT_U(i,:))>0.5))>=2
            mdl = fitlm(GT_U(i,(abs(GT_U(i,:))>0.5)),GT_Y(j+1,(abs(GT_U(i,:))>0.5)));
            R = mdl.Rsquared.Ordinary;
            plot([-4,4],mdl.Coefficients.Estimate(2)*[-4,4]+mdl.Coefficients.Estimate(1),'color','c','LineWidth',0.75);
            annotation('textbox', [0.15, 0.8, 0.35, 0.11],...
            'String', strcat("R-squared: ",num2str(R)), ...
            'FontSize',9);
            annotation('textbox', [0.15, 0.6, 0.23, 0.19],...
            'String', strcat("Coefficients: [",num2str(round(mdl.Coefficients.Estimate(2),4)),",", num2str(round(mdl.Coefficients.Estimate(1),4)),"]"), ...
            'FontSize',9);
        end
       
        plot(GT_U(i,(abs(GT_U(i,:))<=0.5)),GT_Y(j+1,(abs(GT_U(i,:))<=0.5)),'yx','MarkerSize',3); 
        plot([-4,4],[0, 0],'LineWidth',1,'color',[17 17 17]/255);
        plot([0,0],[-1.25, 1.25],'LineWidth',1,'color',[17 17 17]/255);
        f.Position = [620 520 400 200];
        grid on;
        
        hold off;
        title(strcat('Node point','{ }',num2str(i),' - curvature','{ }',num2str(j)));
        xlabel('Normalized curvature (1/m) to 0.001 1/m');
        ylabel('Lateral offset (m)');
        xlim([-4,4]); ylim([-1.25,1.25]);

        xticks([-4:0.5:4]); yticks([-1:0.25:1]);
        
        
        savefig(f, fullfile(temp_folder_path, plots_folder_name,...
                strcat('Correlation_curvature_',num2str((i-1)*3+j),'_',token,'.fig')));
        saveas(f, fullfile(temp_folder_path, plots_folder_name,...
                strcat('Correlation_curvature_',num2str((i-1)*3+j),'_',token,'.png')));
        
        close(f);
    end
end

%% Curvature gradient
if (~isfolder(fullfile(temp_folder_path, plots_folder_name,'curvatureGradient')))
    mkdir(fullfile(temp_folder_path, plots_folder_name,'curvatureGradient'));
end
plots_folder_name = fullfile(plots_folder_name,'curvatureGradient');

% for i=1:size(GT_U,1)
%     for j=1:size(GT_Y,1)-1
%         f = figure();
%         plot(GT_U(i+3,:),GT_Y(j+1,:),'kx','MarkerSize',3);
%         f.Position = [620 520 400 200];
%         hold on;
%         if numel(find(abs(GT_U(i+3,:)) > 0.75)) > 2
%             mdl = fitlm(GT_U(i+3,abs(GT_U(i+3,:)) > 0.75),GT_Y(j+1,abs(GT_U(i+3,:)) > 0.75));
%             R = mdl.Rsquared.Ordinary;
%             plot([-4,0,4],mdl.Coefficients.Estimate(2)*[-4,0,4]+mdl.Coefficients.Estimate(1),'color','b','LineWidth',1.5);
%             annotation('textbox', [0.67, 0.8, 0.35, 0.11],...
%             'String', strcat("R-squared: ",num2str(R)), ...
%             'FontSize',9);
%             annotation('textbox', [0.67, 0.6, 0.23, 0.19],...
%             'String', strcat("Coefficients: [",num2str(round(mdl.Coefficients.Estimate(2),4)),",", num2str(round(mdl.Coefficients.Estimate(1),4)),"]"), ...
%             'FontSize',9);
%         end
%         
%         plot([-4,4],[0, 0],'LineWidth',1,'color',[17 17 17]/255);
%         plot([0,0],[-1.25, 1.25],'LineWidth',1,'color',[17 17 17]/255);
%         
%         grid on;
%         
%         hold off;
%         title(strcat('Node point','{ }',num2str(j),' - curvature gradient', '{ }', num2str(i)));
%         xlabel('Normalized curvature gradient (1/m^2) to 1e-5 1/m^2');
%         ylabel('Lateral offset (m)');
%         xlim([-4,4]); ylim([-1.25,1.25]);
%         xticks([-4:0.5:4]); yticks([-1.25:0.25:1.25]);
%         savefig(f, fullfile(temp_folder_path, plots_folder_name,...
%                 strcat('Correlation_curvaturegradient_',num2str((i-1)*3+j),'_',token,'.fig')));
%         saveas(f, fullfile(temp_folder_path, plots_folder_name,...
%                 strcat('Correlation_curvaturegradient_',num2str((i-1)*3+j),'_',token,'.png')));
%             
%         close(f);
%     end
% end

set(0,'DefaultFigureVisible','on');

end

