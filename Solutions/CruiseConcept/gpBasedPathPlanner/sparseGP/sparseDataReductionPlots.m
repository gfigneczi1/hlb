close all;
clear;

USE_NRMS = true;
if (USE_NRMS)
    ylims = [0, 0.1];
    labelBase = "$NRMS_{e_\delta}";
else
    ylims = [0, 0.6];
    labelBase = "$RMS_{e_\delta}";
end

matFiles = dir(fullfile("C:\git\KDP\publications\GP\results\sparseGP_fixedAcceleration\12drivers\kpis", "KPI_*.mat"));
matFilesEta = dir(fullfile("C:\git\KDP\publications\GP\results\sparseGP_fixedAcceleration\12drivers\kpis", "ETA_*.mat"));

for i=1:length(matFiles)
    data = load(fullfile(matFiles(i).folder, matFiles(i).name));
    driverID = str2num(matFiles(i).name(strfind(matFiles(i).name, 'driver_')+7:strfind(matFiles(i).name, '.')-1));
    inputID = str2num(matFiles(i).name(strfind(matFiles(i).name, 'indSize_')+8:strfind(matFiles(i).name, '_driver_')-1));
    if (inputID > 0)
        dataEta = load(fullfile(matFilesEta(i).folder, matFilesEta(i).name));
        N = size(dataEta.ETA(1).input_estimation, 2);
        M = numel(dataEta.ETA(1).hyp_opt.cov);
        for shiftID = 1:numel(data.KPI)
            if(USE_NRMS)
                NRMS_val(driverID, inputID, shiftID) = data.KPI{shiftID}(2);
                NRMS_tr(driverID, inputID, shiftID) = data.KPI{shiftID}(6);
                NRMS_ind(driverID, inputID, shiftID) = data.KPI{shiftID}(8+M+2);
                timeRun(driverID, inputID, shiftID) = data.KPI{shiftID}(8+M+4);
            else
                % using RMS upscaling by the central normalization factor
                eta = load(fullfile(matFilesEta(i).folder, matFilesEta(i).name));
                sigma = eta.ETA(shiftID).normFactors.s_out(shiftID);
    
                NRMS_val(driverID, inputID, shiftID) = data.KPI{shiftID}(1)*sigma;
                NRMS_tr(driverID, inputID, shiftID) = data.KPI{shiftID}(5)*sigma;
                NRMS_ind(driverID, inputID, shiftID) = data.KPI{shiftID}(8+M+1)*sigma;
                timeRun(driverID, inputID, shiftID) = data.KPI{shiftID}(8+M+4);
            end
        end
    else
        for shiftID = 1:numel(data.KPI)
            NRMS_val_ref(driverID, 1, shiftID) = data.KPI{shiftID}(2);
            NRMS_tr_ref(driverID, 1, shiftID) = data.KPI{shiftID}(6);
            NRMS_ind_ref(driverID, 1, shiftID) = data.KPI{shiftID}(27);
            timeRun_ref(driverID, 1, shiftID) = data.KPI{shiftID}(29);
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

NRMS_ind_avgd_ref = [];
NRMS_tr_avgd_ref = [];
NRMS_ind_avgd_ref = [];
timeRun_avgd_ref = [];

for driverID=1:size(NRMS_val_ref,1)
    for inputID=1:size(NRMS_val_ref,2)
        NRMS_val_avgd_ref(driverID,inputID) = mean(NRMS_val_ref(driverID,inputID,:));
        NRMS_tr_avgd_ref(driverID,inputID) = mean(NRMS_tr_ref(driverID,inputID,:));
        NRMS_ind_avgd_ref(driverID,inputID) = mean(NRMS_ind_ref(driverID,inputID,:));
        timeRun_avgd_ref(driverID,inputID) = mean(timeRun_ref(driverID,inputID,:));
    end
end

markers = {'x', 'o', '*', '.', '+', 'diamond', 'square', '^', '<', '>', 'hexagram', 'pentagram', 'none', 'x'};
variables = ["$10$", "$20$", "$50$", "$100$", "$300$", "$500$", "$1000$", "FullGP"];

%% DRIVERS' data plotted in the function of the different epsilon values - induction, train and validation data
f = figure(1);
f.Position = [100 100 1050 750];
set(f,'defaulttextInterpreter','latex') ;
set(f, 'defaultAxesTickLabelInterpreter','latex');  
set(f, 'defaultLegendInterpreter','latex');

subplot(3,3,1); % induction plots

for driverID=1:size(NRMS_ind_avgd,1)
    plot([NRMS_ind_avgd(driverID,:) NRMS_ind_avgd_ref(driverID,:)], 'Marker', markers{driverID}, 'DisplayName', strcat("Driver", {' '}, num2str(driverID)));
    grid on;
    hold on;
end
%legend("Location", "southoutside", "NumColumns", 5);
set(gca,'FontSize', 14);
title("Induction accuracy");
ylabel(strcat(labelBase,"^{ind}$"));
xlabel("$\varepsilon$");
ylim(ylims);


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
xlabel("$\varepsilon$");

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
xlabel("$\varepsilon$");
ylim(ylims);

xticks([1:1:size(NRMS_tr_avgd,2)+1]);
xticklabels([]);

set(gca,'TickLabelInterpreter','latex');

subplot(3,3,5);
boxplot([NRMS_tr_avgd(:,1:end) NRMS_tr_avgd_ref(:,1:end)]);
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
    plot([NRMS_val_avgd(driverID,:) NRMS_val_avgd_ref(driverID,:)], 'Marker', markers{driverID}, 'DisplayName', strcat("Driver", {' '}, num2str(driverID)));
    grid on;
    hold on;
end
%legend("Location", "southoutside", "NumColumns", 5);
set(gca,'FontSize', 14);
title("Validation accuracy");
ylabel(strcat(labelBase,"^{val}$"));
xlabel("$\varepsilon$");
ylim(ylims);

xticks([1:1:size(NRMS_val_avgd,2)+1]);
xticklabels([]);

set(gca,'TickLabelInterpreter','latex');

subplot(3,3,6);
boxplot([NRMS_val_avgd(:,1:end) NRMS_val_avgd_ref(:,1:end)]);
xticklabels([]);
grid on;
ylim(ylims);

set(gca,'FontSize', 14);
title("Validation accuracy");
ylabel(strcat(labelBase,"^{val}$"));
xlabel("$\varepsilon$");

set(gca,'TickLabelInterpreter','latex');

subplot(3,1,3); % runtime plots
boxplot([timeRun_avgd(:,1:end) timeRun_avgd_ref(:,1:end)]);
xticklabels(variables(1:end));
grid on;

set(gca,'FontSize', 14);
title("Run time");
ylabel("$T_c(s)$");
xlabel("$\varepsilon$");

xticks([1:1:size(timeRun_avgd,2)+1]);

set(gca,'TickLabelInterpreter','latex');

