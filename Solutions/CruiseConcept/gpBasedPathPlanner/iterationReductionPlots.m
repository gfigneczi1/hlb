close all;
clear;

USE_NRMS = false;
TYPE = "PP3";

if (USE_NRMS)
    ylims = [0, 0.15];
    labelBase = "$NRMS_{e_\delta}";
else
    ylims = [0, 0.15];
    labelBase = "$RMS_{e_\delta}";
end


matFiles = dir(fullfile("C:\git\KDP\publications\GP\results\iteration_reduction_withPP3\KPI_input*.mat"));

for i=1:length(matFiles)
    data = load(fullfile(matFiles(i).folder, matFiles(i).name));
    driverID = str2num(matFiles(i).name(strfind(matFiles(i).name, 'driver_')+7:strfind(matFiles(i).name, '.')-1));
    inputID = str2num(matFiles(i).name(strfind(matFiles(i).name, 'input_')+6:strfind(matFiles(i).name, '_driver_')-1));
    for shiftID = 1:numel(data.KPI)
        if(USE_NRMS)
            NRMS_val(driverID, inputID, shiftID) = data.KPI{shiftID}(2);
            NRMS_tr(driverID, inputID, shiftID) = data.KPI{shiftID}(6);
            NRMS_ind(driverID, inputID, shiftID) = data.KPI{shiftID}(27);
            timeRun(driverID, inputID, shiftID) = data.KPI{shiftID}(29);
        else
            % using RMS upscaling by the central normalization factor
            etaFiles = dir(fullfile("C:\git\KDP\publications\GP\results\iteration_reduction_withPP3\ETA_input*.mat"));
            eta = load(fullfile(etaFiles(i).folder, etaFiles(i).name));
            sigma = eta.ETA(shiftID).normFactors(17+shiftID);

            NRMS_val(driverID, inputID, shiftID) = data.KPI{shiftID}(1)*sigma;
            NRMS_tr(driverID, inputID, shiftID) = data.KPI{shiftID}(5)*sigma;
            NRMS_ind(driverID, inputID, shiftID) = data.KPI{shiftID}(26)*sigma;
            timeRun(driverID, inputID, shiftID) = data.KPI{shiftID}(29);
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
variables = ["$5$", "$10$", "$20$", "$50$", "$100$"];

%% DRIVERS' data plotted in the function of the different epsilon values - induction, train and validation data
f = figure(1);
f.Position = [100 100 1050 750];
set(f,'defaulttextInterpreter','latex') ;
set(f, 'defaultAxesTickLabelInterpreter','latex');  
set(f, 'defaultLegendInterpreter','latex');

subplot(3,3,1); % induction plots
for driverID=1:size(NRMS_ind_avgd,1)
    plot([NRMS_ind_avgd(driverID,:) NRMS_ind_avgd(driver, 'Marker', markers{driverID}, 'DisplayName', strcat("Driver", {' '}, num2str(driverID)));
    grid on;
    hold on;
end
%legend("Location", "southoutside", "NumColumns", 5);
set(gca,'FontSize', 14);
title("Induction accuracy");
ylabel(strcat(labelBase,"^{ind}$"));
xlabel("$N$");

xticks([1:1:size(NRMS_ind_avgd,2)+1]);
xticklabels([]);

set(gca,'TickLabelInterpreter','latex');

subplot(3,3,4);
boxplot([NRMS_ind_avgd(:,1:end) NRMS_ind_avgd_ref(:,1:end)]);
xticklabels([]);
grid on;
ylim(ylims);

set(gca,'FontSize', 14);
title("Induction accuracy");
ylabel(strcat(labelBase,"^{ind}$"));
xlabel("$N$");

set(gca,'TickLabelInterpreter','latex');

subplot(3,3,2); % train data accuracy
for driverID=1:size(NRMS_tr_avgd,1)
    plot([NRMS_tr_avgd(driverID,:) NRMS_tr_avgd_ref(driverID,:)], 'Marker', markers{driverID}, 'DisplayName', strcat("Driver", {' '}, num2str(driverID)));
    grid on;
    hold on;
end
%legend("Location", "southoutside", "NumColumns", 5);
set(gca,'FontSize', 14);
title("Train accuracy");
ylabel(strcat(labelBase,"^{tr}$"));
xlabel("$N$");
ylim(ylims);

xticks([1:1:size(NRMS_tr_avgd,2)+1]);
xticklabels([]);

set(gca,'TickLabelInterpreter','latex');

subplot(3,3,5);
boxplot([NRMS_tr_avgd(:,1:end) NRMS_tr_avgd_ref]);
xticklabels([]);
grid on;
ylim(ylims);

set(gca,'FontSize', 14);
title("Train accuracy");
ylabel(strcat(labelBase,"^{tr}$"));
xlabel("$N$");

set(gca,'TickLabelInterpreter','latex');

subplot(3,3,3); % validation data accuracy
for driverID=1:size(NRMS_val_avgd,1)
    plot([NRMS_val_avgd(driverID,:) NRMS_val_avgd_ref(driverID,:)], 'Marker', markers{driverID}, 'DisplayName', strcat("Driver", {' '}, num2str(driverID)));
    grid on;
    hold on;
end
%legend("Location", "southoutside", "NumColumns", 5);
set(gca,'FontSize', 14);
title("Validation accuracy");
ylabel(strcat(labelBase,"^{val}$"));
xlabel("$N$");
ylim(ylims);

xticks([1:1:size(NRMS_val_avgd,2)+1]);
xticklabels([]);

set(gca,'TickLabelInterpreter','latex');

subplot(3,3,6);
boxplot([NRMS_val_avgd(:,1:end) NRMS_val_avgd_ref]);
xticklabels([]);
grid on;
ylim(ylims);

set(gca,'FontSize', 14);
title("Validation accuracy");
ylabel(strcat(labelBase,"^{val}$"));
xlabel("$N$");

set(gca,'TickLabelInterpreter','latex');

subplot(3,1,3); % runtime plots
boxplot([timeRun_avgd(:,1:end) timeRun_avgd_ref]);
xticklabels(variables(1:end));
grid on;

set(gca,'FontSize', 14);
title("Run time");
ylabel("$T_c(s)$");
xlabel("$N$");

xticks([1:1:size([timeRun_avgd timeRun_avgd_ref],2)]);

set(gca,'TickLabelInterpreter','latex');

