function plot_driverCorrelationWithModel(config,GT_U_Driver, GT_Y_Driver, GT_U_model, GT_Y_model, token)

global model_ver

temp_folder_path = config.root;
plots_folder_name = 'plots';
set(0,'DefaultFigureVisible','off');

if (model_ver == 1)
    for i=1:size(GT_U_Driver,1)-4
        for j=1:size(GT_Y_Driver,1)-1
            f = figure();
            plot(GT_U_Driver(i,(GT_U_Driver(i,:)>0)),GT_Y_Driver(j+1,(GT_U_Driver(i,:)>0)),'kx','MarkerSize',3); 
            hold on;
            plot(GT_U_model(i,(GT_U_model(i,:))>0),GT_Y_model(j+1,(GT_U_model(i,:)>0)),'bo','MarkerSize',3); 
            legend('Driver', 'Model');

            plot(GT_U_Driver(i+3,(GT_U_Driver(i+3,:))<0),GT_Y_Driver(j+1,(GT_U_Driver(i+3,:)<0)),'kx','MarkerSize',3); 
            plot(GT_U_model(i+3,(GT_U_model(i+3,:)<0)),GT_Y_model(j+1,(GT_U_model(i+3,:)<0)),'bo','MarkerSize',3); 

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
                    strcat('Correlation_comparison_',num2str((i-1)*3+j),'_',token,'.fig')));
            saveas(f, fullfile(temp_folder_path, plots_folder_name,...
                    strcat('Correlation_comparison_',num2str((i-1)*3+j),'_',token,'.png')));

            close(f);
        end
    end
end

set(0,'DefaultFigureVisible','on');

end

