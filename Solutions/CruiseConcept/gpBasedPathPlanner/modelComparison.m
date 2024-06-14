close all;
clear;

matFiles = dir(fullfile("C:\git\KDP\publications\GP\results\linearRegression_model\KPI_LDM*.mat"));
mode = "train"; %validation

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
        NRMS_LDM(driverID, shiftID) = data.KPI(shiftID,idLDM);
    end
    NRMS_avg_LDM(driverID) = mean(NRMS_LDM(driverID, :));
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

plot([1:1:length(NRMS_avg_GP)], NRMS_avg_GP, 'Marker', 'x', 'color', 'g');
hold on; grid on;
plot([1:1:length(NRMS_avg_LDM)], NRMS_avg_LDM, 'Marker', 'o', 'color', 'b');
if (mode=="train")
    legend("$NRMS_{e_{\delta},GPM}^{tr}$", "$NRMS_{e_{\delta},LRM}^{tr}$");
    title("Model estimation comparison - train data");
else
    legend("$NRMS_{e_{\delta},GPM}^{val}$", "$NRMS_{e_{\delta},LRM}^{val}$", 'Location', 'best');
    title("Model estimation comparison - validation data");
end
ylabel("$$NRMS_{e_{\delta}}$"); xlabel("Driver ID");
xticks([1:1:length(NRMS_GP)]);
ylim([0, 1.1*max(max(NRMS_avg_GP),max(NRMS_avg_LDM))]);

set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize', 14);