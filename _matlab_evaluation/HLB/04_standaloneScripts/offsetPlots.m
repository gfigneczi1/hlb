clear all
close all
sourceOfInformation = "Resim";
%folderLDM = 'C:\git\KDP_Igneczi\publikációk\ELDM\plots_ELDM_overfitted\plots';
folderLDM = 'C:\Users\igg2bp\Documents\00_Projects\Research\LDM_comparison_article';
if (sourceOfInformation == "Measurement")

    files = dir(fullfile(folderLDM,'Dr*Model*.mat'));

    for i=1:length(files)
        load(fullfile(files(i).folder, files(i).name));
        drID = str2num(files(i).name(3:5));
        curveID = files(i).name(end-5:end-4);
        if (contains(curveID,'_'))
            curveID = curveID(2);
        end
        curveID = str2num(curveID);
        if (contains(files(i).name,"62B"))
            d{drID, 10-curveID+1} = KPI.curves.data;
        else
            d{drID, curveID} = KPI.curves.data;
        end
    end
elseif (sourceOfInformation == "Resim")
    d = load('C:\Users\igg2bp\Documents\00_Projects\Research\LDM_comparison_article\data.mat');
    d = d.data;
end

% 1 dimensional classification - only drivers are classified
% based on driver ID
drClasses = [1 2 2 2 3 ...
    1 2 3 1 2 ...
    3 3 2 2 1];
classIDs = unique(drClasses);
for j=1:size(d,2)
    for k=1:length(classIDs)
        f = figure(); 
        for i=1:size(d,1)
            if (drClasses(i)==classIDs(k))
                if(~isempty(d{i,j}))
                    traj_y = spline(d{i,j}.traj(:,1), d{i,j}.traj(:,2), d{i,j}.cor(:,1));
                    ref_y = spline(d{i,j}.ref(:,1), d{i,j}.ref(:,2), d{i,j}.cor(:,1));
                    % global to local transformation to calculate the
                    % offset
                    traj = [d{i,j}.cor(:,1) traj_y];
                    cor = [d{i,j}.cor(:,1) d{i,j}.cor(:,2)];
                    ref = [d{i,j}.ref(:,1) ref_y];
                    T = [cos(d{i,j}.orient(1)) sin(d{i,j}.orient(1)); -sin(d{i,j}.orient(1)) cos(d{i,j}.orient(1))];
                    traj_local = traj;
                    cor_local = cor;
                    ref_local = ref;
                    offsetsTraj = traj_local(:,2) - cor_local(:,2);                    
                    offsetsRef = ref_local(:,2) - cor_local(:,2);
                    if (((d{i,j}.ref(end,1)-d{i,j}.ref(1,1)).^2+(d{i,j}.ref(end,2)-d{i,j}.ref(1,2)).^2).^0.5 > 1e4)
                        measurementUncut = 1;
                    else
                        measurementUncut = 0;
                    end
                    subplot(4,1,1);
                    hold on;
                    title(strcat('Driver class', {' '}, num2str(k), ' curve', {' '}, num2str(j)));
                    
                    if (~measurementUncut)
                        plot(d{i,j}.cor(:,1), offsetsTraj, 'DisplayName', strcat('Dr',num2str(i)));
                    else
                        plot([],[],'DisplayName', strcat('Dr',num2str(i), '-missing'));
                    end
                     ylim([-1, 1]); 
                    grid on;
                    xlabel('X(m)'); ylabel('offset(m)');
                    legend;
                    subplot(4,1,2);
                    hold on;
                    if (~measurementUncut)
                        plot(d{i,j}.cor(:,1), offsetsRef,  'DisplayName', strcat('Dr',num2str(i)));
                    else
                        plot([],[],'DisplayName', strcat('Dr',num2str(i), '-missing'));
                    end
                    grid on;
                    xlabel('X(m)'); ylabel('offset(m)');
                    ylim([-1, 1]); 
                    legend;
                    subplot(4,1,3);
                    plot(d{i,j}.ref(:,1), d{i,j}.ref(:,2));
                    grid on;
                    xlabel('X(m)');
                    ylabel('Y(m)');
                    subplot(4,1,4);
                    plot(d{i,j}.ref(:,1), d{i,j}.curv(:,1));
                    grid on;
                    title('Curvature');
                    xlabel('X(m)'); ylabel('curvature(1/m)');
                end
            end
        end
        saveas(f, fullfile(folderLDM,...
                        strcat('Offsets_DriverClass_',num2str(classIDs(k)), '_Curve_', num2str(j),'.png')));
        savefig(f, fullfile(folderLDM,...
              strcat('Offsets_DriverClass_',num2str(classIDs(k)), '_Curve_', num2str(j),'.fig')));
    	close(f);
    end
end

