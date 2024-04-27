close all;
clear;

matFiles = dir(fullfile("C:\git\KDP\publications\GP\results\input_reduction\KPI_input*.mat"));
id = 6;

for i=1:length(matFiles)
    data = load(fullfile(matFiles(i).folder, matFiles(i).name));
    driverID = str2num(matFiles(i).name(strfind(matFiles(i).name, 'driver_')+7:strfind(matFiles(i).name, '.')-1));
    inputID = str2num(matFiles(i).name(strfind(matFiles(i).name, 'input_')+6:strfind(matFiles(i).name, '_driver_')-1))+1;
    for shiftID = 1:numel(data.KPI)
        NRMS(driverID, inputID, shiftID) = data.KPI{shiftID}(id);
    end
end

%% NOW DOING PLOTS FOR ALL INPUTS SEPARATELY

for driverID=1:size(NRMS,1)
    for inputID=1:size(NRMS,2)
        NRMS_avgd(driverID,inputID) = mean(NRMS(driverID,inputID,:));
    end
    NRMS2D{driverID}(:,:) = NRMS(driverID,:,:);
end

markers = {'x', 'o', '*', '.', '+', 'diamond', 'square', '^', '<', '>', 'hexagram', 'pentagram', 'none', 'x'};
variables = ["$t_{pass}$", "$o_{t}$", "$fO_{t}$", "$v_x$", "$a_x$", "$\omega$", "$\kappa_{road}$", "$d\kappa$", "$t_{pass}$ and $o_{t}$","$\omega$ and $\kappa$", "$t_{pass}$, $fO_t$ and $o_{t}$", "$fO_t$ and $O_t$", "$fO_t$ and $\omega_z$"];

f = figure(1);
f.Position = [100 100 850 500];
set(f,'defaulttextInterpreter','latex') ;
set(f, 'defaultAxesTickLabelInterpreter','latex');  
set(f, 'defaultLegendInterpreter','latex');

inputID = 1;
plot(NRMS_avgd(:,inputID), 'Marker', markers{inputID}, 'LineWidth', 2, 'DisplayName', "All inputs");
    grid on;
    hold on;
for inputID=10:size(NRMS_avgd,2)
    if (inputID==1)
        plot(NRMS_avgd(:,inputID), 'Marker', markers{inputID}, 'LineWidth', 2, 'DisplayName', "All inputs");
    else
        plot(NRMS_avgd(:,inputID), 'Marker', markers{inputID}, 'DisplayName', strcat("No input", {' '}, variables(inputID-1)));
    end
    grid on;
    hold on;
    xlabel("Driver ID");
end

legend("Location", "southoutside", "NumColumns", 3);
set(gca,'FontSize', 14);
if (id==2)
    title("Estimation accuracy degredation - validation data");
    ylabel("$NRMS_{e_\delta}^{val}$");
else
    title("Estimation accuracy degredation - train data");
    ylabel("$NRMS_{e_\delta}^{tr}$");
end
ylim([0, max(max(NRMS_avgd))*1.1]);

f = figure(2);
set(f,'defaulttextInterpreter','latex') ;
set(f, 'defaultAxesTickLabelInterpreter','latex');  
set(f, 'defaultLegendInterpreter','latex');

for driverID=1:9 %numel(NRMS2D)
    subplot(2,5,driverID);
    surf(NRMS2D{driverID}');
    xlabel("inputs");
    ylabel("nodePoint");
    xlim([1,14]);
    xticks(1:1:14);
    xticklabels(variables);
    title(strcat("Driver", num2str(driverID)));
    ylim([1,10]);
    view(2);
end
