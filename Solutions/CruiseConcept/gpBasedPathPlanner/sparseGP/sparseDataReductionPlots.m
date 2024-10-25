close all;
clear;

USE_NRMS = true;
if (USE_NRMS)
    ylims = [0, 0.25];
    labelBase = "$NRMS_{e_\delta}";
else
    ylims = [0, 0.6];
    labelBase = "$RMS_{e_\delta}";
end

matFiles = dir(fullfile("kpis", "KPI_indSize*.mat"));

for i=1:length(matFiles)
    data = load(fullfile(matFiles(i).folder, matFiles(i).name));
    driverID = str2num(matFiles(i).name(strfind(matFiles(i).name, 'driver_')+7:strfind(matFiles(i).name, '.')-1));
    inputID = str2num(matFiles(i).name(strfind(matFiles(i).name, 'indSize_')+8:strfind(matFiles(i).name, '_driver_')-1));
    for shiftID = 1:numel(data.KPI)
        if(USE_NRMS)
            NRMS_val(driverID, inputID, shiftID) = data.KPI{shiftID}(2);
            NRMS_tr(driverID, inputID, shiftID) = data.KPI{shiftID}(6);
            NRMS_ind(driverID, inputID, shiftID) = data.KPI{shiftID}(18);
            timeRun(driverID, inputID, shiftID) = data.KPI{shiftID}(20);
        else
            % using RMS upscaling by the central normalization factor
            etaFiles = dir(fullfile("kpis", "ETA_indSize*.mat"));
            eta = load(fullfile(etaFiles(i).folder, etaFiles(i).name));
            sigma = eta.ETA.normFactors(10+shiftID);

            NRMS_val(driverID, inputID, shiftID) = data.KPI{shiftID}(1)*sigma;
            NRMS_tr(driverID, inputID, shiftID) = data.KPI{shiftID}(5)*sigma;
            NRMS_ind(driverID, inputID, shiftID) = data.KPI{shiftID}(17)*sigma;
            timeRun(driverID, inputID, shiftID) = data.KPI{shiftID}(20);
        end
    end
end

%% NOW DOING PLOTS FOR ALL INPUTS SEPARATELY

for driverID=1:size(NRMS_val,1)
    for inputID=1:size(NRMS_val,2)
        NRMS_val_avgd(driverID,inputID) = mean(NRMS_val(driverID,inputID,:));
        NRMS_tr_avgd(driverID,inputID) = mean(NRMS_tr(driverID,inputID,:));
        NRMS_ind_avgd(driverID,inputID) = mean(NRMS_ind(driverID,inputID,:));
        timeRun_avgd(driverID,inputID) = mean(timeRun(driverID,inputID,:));
    end
end

markers = {'x', 'o', '*', '.', '+', 'diamond', 'square', '^', '<', '>', 'hexagram', 'pentagram', 'none', 'x'};
variables = ["$10$", "$20$", "$50$", "$100$", "$300$"];

%% DRIVERS' data plotted in the function of the different epsilon values - induction, train and validation data
f = figure(1);
f.Position = [100 100 1050 750];
set(f,'defaulttextInterpreter','latex') ;
set(f, 'defaultAxesTickLabelInterpreter','latex');  
set(f, 'defaultLegendInterpreter','latex');

subplot(3,3,1); % induction plots
for driverID=1:size(NRMS_ind_avgd,1)
    plot(NRMS_ind_avgd(driverID,:), 'Marker', markers{driverID}, 'DisplayName', strcat("Driver", {' '}, num2str(driverID)));
    grid on;
    hold on;
end
%legend("Location", "southoutside", "NumColumns", 5);
set(gca,'FontSize', 14);
title("Induction accuracy");
ylabel(strcat(labelBase,"^{ind}$"));
xlabel("$\varepsilon$");
ylim(ylims);


xticks([1:1:size(NRMS_ind_avgd,2)]);
xticklabels([]);

set(gca,'TickLabelInterpreter','latex');

subplot(3,3,4);
boxplot(NRMS_ind_avgd(:,1:end));
xticklabels([]);
grid on;
ylim(ylims);

set(gca,'FontSize', 14);
title("Induction accuracy");
ylabel(strcat(labelBase,"^{ind}$"));
xlabel("$\varepsilon$");

set(gca,'TickLabelInterpreter','latex');

subplot(3,3,2); % train data accuracy
for driverID=1:size(NRMS_tr_avgd,1)
    plot(NRMS_tr_avgd(driverID,:), 'Marker', markers{driverID}, 'DisplayName', strcat("Driver", {' '}, num2str(driverID)));
    grid on;
    hold on;
end
%legend("Location", "southoutside", "NumColumns", 5);
set(gca,'FontSize', 14);
title("Train accuracy");
ylabel(strcat(labelBase,"^{tr}$"));
xlabel("$\varepsilon$");
ylim(ylims);

xticks([1:1:size(NRMS_tr_avgd,2)]);
xticklabels([]);

set(gca,'TickLabelInterpreter','latex');

subplot(3,3,5);
boxplot(NRMS_tr_avgd(:,1:end));
xticklabels([]);
grid on;
ylim(ylims);

set(gca,'FontSize', 14);
title("Train accuracy");
ylabel(strcat(labelBase,"^{tr}$"));
xlabel("$\varepsilon$");

set(gca,'TickLabelInterpreter','latex');

subplot(3,3,3); % validation data accuracy
for driverID=1:size(NRMS_val_avgd,1)
    plot(NRMS_val_avgd(driverID,:), 'Marker', markers{driverID}, 'DisplayName', strcat("Driver", {' '}, num2str(driverID)));
    grid on;
    hold on;
end
%legend("Location", "southoutside", "NumColumns", 5);
set(gca,'FontSize', 14);
title("Validation accuracy");
ylabel(strcat(labelBase,"^{val}$"));
xlabel("$\varepsilon$");
ylim(ylims);

xticks([1:1:size(NRMS_val_avgd,2)]);
xticklabels([]);

set(gca,'TickLabelInterpreter','latex');

subplot(3,3,6);
boxplot(NRMS_val_avgd(:,1:end));
xticklabels([]);
grid on;
ylim(ylims);

set(gca,'FontSize', 14);
title("Validation accuracy");
ylabel(strcat(labelBase,"^{val}$"));
xlabel("$\varepsilon$");

set(gca,'TickLabelInterpreter','latex');

subplot(3,1,3); % runtime plots
boxplot(timeRun_avgd(:,1:end));
xticklabels(variables(1:end));
grid on;

set(gca,'FontSize', 14);
title("Run time");
ylabel("$T_c(s)$");
xlabel("$\varepsilon$");

xticks([1:1:size(timeRun_avgd,2)]);

set(gca,'TickLabelInterpreter','latex');

