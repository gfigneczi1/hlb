close all;
clear;

matFiles = dir(fullfile("C:\git\KDP\publications\GP\results\input_reduction\KPI_input*.mat"));
id = 2; %6: train data, 2: validation data

for i=1:length(matFiles)
    data = load(fullfile(matFiles(i).folder, matFiles(i).name));
    driverID = str2num(matFiles(i).name(strfind(matFiles(i).name, 'driver_')+7:strfind(matFiles(i).name, '.')-1));
    inputID = str2num(matFiles(i).name(strfind(matFiles(i).name, 'input_')+6:strfind(matFiles(i).name, '_driver_')-1))+1;
    for shiftID = 1:numel(data.KPI)
        NRMS(driverID, inputID, shiftID) = data.KPI{shiftID}(id);
    end
end

%% NOW DOING PLOTS FOR ALL INPUTS SEPARATELY

for driverID=1:8
    for inputID=1:size(NRMS,2)
        NRMS_avgd(driverID,inputID) = mean(NRMS(driverID,inputID,:));
    end
    NRMS2D{driverID}(:,:) = NRMS(driverID,:,:);
end

markers = {'x', 'o', '*', '.', '+', 'diamond', 'square', '^', '<', '>', 'hexagram', 'pentagram', 'none', 'x'};
variables = ["$t_{pass}$", "$o_{t}$", "$fO_{t}$", "$v_x$", "$a_x$", "$\omega$", "$\kappa_{road}$", "$d\kappa$", "$fO_t$ and $o_{t}$", "$t_{pass}$, $fO_t$ and $o_{t}$", "$t_{pass}$ and $fO_t$", "$fO_t$ and $d\kappa$"];

f = figure(1);
f.Position = [100 100 850 500];
set(f,'defaulttextInterpreter','latex') ;
set(f, 'defaultAxesTickLabelInterpreter','latex');  
set(f, 'defaultLegendInterpreter','latex');

inputID = 1;
plot(NRMS_avgd(:,inputID), 'Marker', markers{inputID}, 'LineWidth', 2, 'DisplayName', "All inputs");
    grid on;
    hold on;
for inputID=2:size(NRMS_avgd,2)
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

%% DRIVERS' data plotted in the function of the inputs
f = figure(4);
f.Position = [100 100 850 650];
set(f,'defaulttextInterpreter','latex') ;
set(f, 'defaultAxesTickLabelInterpreter','latex');  
set(f, 'defaultLegendInterpreter','latex');

subplot(2,1,1);
for driverID=1:size(NRMS_avgd,1)
    plot(NRMS_avgd(driverID,:), 'Marker', markers{driverID}, 'DisplayName', strcat("Driver", {' '}, num2str(driverID)));
    grid on;
    hold on;
end
legend("Location", "southoutside", "NumColumns", 5);
set(gca,'FontSize', 14);
if (id==2)
    title("Estimation accuracy degredation - validation data");
    ylabel("$NRMS_{e_\delta}^{val}$");
else
    title("Estimation accuracy degredation - train data");
    ylabel("$NRMS_{e_\delta}^{tr}$");
end
ylim([0, max(max(NRMS_avgd))*1.1]);
xticks([1:1:size(NRMS_avgd,2)]);
xticklabels(["All inputs", variables(1:size(NRMS_avgd,2)-1)]);

subplot(2,1,2);
boxplot(NRMS_avgd(:,1:end));
xticklabels(["All inputs", variables(1:end)]);
grid on;
ylim([0,0.15]);

set(gca,'FontSize', 14);
if (id==2)
    title("Estimation accuracy degredation - validation data");
    ylabel("$NRMS_{e_\delta}^{val}$");
else
    title("Estimation accuracy degredation - train data");
    ylabel("$NRMS_{e_\delta}^{tr}$");
end

set(gca,'TickLabelInterpreter','latex');

f = figure(2);
set(f,'defaulttextInterpreter','latex') ;
set(f, 'defaultAxesTickLabelInterpreter','latex');  
set(f, 'defaultLegendInterpreter','latex');

for driverID=1:numel(NRMS2D)
    subplot(2,5,driverID);
    surf(NRMS2D{driverID}');
    xlabel("inputs");
    ylabel("nodePoint");
    xlim([1,size(NRMS,2)]);
    xticks(1:1:size(NRMS,2));
    xticklabels(["All input", variables(1:size(NRMS,2)-1)]);
    title(strcat("Driver", num2str(driverID)));
    ylim([1,10]);
    view(2);
end

f = figure(3);
f.Position = [100 100 850 300];
set(f,'defaulttextInterpreter','latex') ;
set(f, 'defaultAxesTickLabelInterpreter','latex');  
set(f, 'defaultLegendInterpreter','latex');

boxplot(NRMS_avgd(:,1:end));
xticklabels(["All inputs", variables(1:end)]);
grid on;
ylim([0,0.15]);

set(gca,'FontSize', 14);
if (id==2)
    title("Estimation accuracy degredation - validation data");
    ylabel("$NRMS_{e_\delta}^{val}$");
else
    title("Estimation accuracy degredation - train data");
    ylabel("$NRMS_{e_\delta}^{tr}$");
end

set(gca,'TickLabelInterpreter','latex');
