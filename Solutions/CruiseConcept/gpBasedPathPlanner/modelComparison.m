close all;
clear;

matFiles = dir(fullfile("C:\git\KDP\publications\GP\results\linearRegression_model\KPI_LDM*.mat"));
mode = "validation"; %validation

if (mode=="train")
    idLDM = 5; %5: train data, 2: validation data
    idGP = 6;
else
    idLDM = 2; idGP = 2;
end

% Linear regression model
for i=1:length(matFiles)
    data = load(fullfile(matFiles(i).folder, matFiles(i).name));
    driverID = str2num(matFiles(i).name(strfind(matFiles(i).name, 'driver_')+7:strfind(matFiles(i).name, '.')-1));
    for shiftID = 1:size(data.KPI,1)
        NRMS_LRM(driverID, shiftID) = data.KPI(shiftID,idLDM);
    end
    NRMS_avg_LRM(driverID) = mean(NRMS_LRM(driverID, :));
end

% ELDM
matFiles = dir("C:\git\KDP\publications\GP\results\eldm\KPI_ELDM*.mat");
for i=1:length(matFiles)
    data = load(fullfile(matFiles(i).folder, matFiles(i).name));
    driverID = str2num(matFiles(i).name(strfind(matFiles(i).name, 'driver_')+7:strfind(matFiles(i).name, '.')-1));
    for shiftID = 1:size(data.KPI,1)
        NRMS_ELDM(driverID, shiftID) = data.KPI(shiftID,idLDM);
    end
    NRMS_avg_ELDM(driverID) = mean(NRMS_ELDM(driverID, :));
end

% GP model
matFiles = dir(fullfile("C:\git\KDP\publications\GP\results\gp_model\KPI_input*.mat"));
for i=1:length(matFiles)
    data = load(fullfile(matFiles(i).folder, matFiles(i).name));
    driverID = str2num(matFiles(i).name(strfind(matFiles(i).name, 'driver_')+7:strfind(matFiles(i).name, '.')-1));
    for shiftID = 1:numel(data.KPI)
        NRMS_GP(driverID, shiftID) = data.KPI{shiftID}(idGP);
    end
    NRMS_avg_GP(driverID) = mean(NRMS_GP(driverID, :));
end

%% create comparison plots
f = figure(1);
f.Position = [100 100 550 350];
set(f,'defaulttextInterpreter','latex') ;
set(f, 'defaultAxesTickLabelInterpreter','latex');  
set(f, 'defaultLegendInterpreter','latex');

plot([1:1:length(NRMS_avg_GP)], NRMS_avg_GP, 'Marker', 'x', 'color', 'r');
hold on; grid on;
plot([1:1:length(NRMS_avg_LRM)], NRMS_avg_LRM, 'Marker', 'o', 'color', 'b');
plot([1:1:length(NRMS_avg_ELDM)], NRMS_avg_ELDM, 'Marker', 'o', 'color', 'g');
if (mode=="train")
    legend("$NRMS_{e_{\delta},GPM}^{tr}$", "$NRMS_{e_{\delta},LRM}^{tr}$", "$NRMS_{e_{\delta},ELDM}^{tr}$");
    title("Model estimation comparison - train data");
else
    legend("$NRMS_{e_{\delta},GPM}^{val}$", "$NRMS_{e_{\delta},LRM}^{val}$", "$NRMS_{e_{\delta},ELDM}^{val}$", 'Location', 'best');
    title("Model estimation comparison - validation data");
end
ylabel("$$NRMS_{e_{\delta}}$"); xlabel("Driver ID");
xticks([1:1:length(NRMS_GP)]);
ylim([0, 1.1*max([max(NRMS_avg_GP),max(NRMS_avg_LRM), max(NRMS_avg_ELDM)])]);

set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize', 14);