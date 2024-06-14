close all;
clear;

matFiles = dir(fullfile("C:\git\KDP\publications\GP\results\linearRegression_model\KPI_LDM*.mat"));
id = 5; %5: train data, 2: validation data

for i=1:length(matFiles)
    data = load(fullfile(matFiles(i).folder, matFiles(i).name));
    driverID = str2num(matFiles(i).name(strfind(matFiles(i).name, 'driver_')+7:strfind(matFiles(i).name, '.')-1));
    for shiftID = 1:size(data.KPI,1)
        NRMS(driverID, shiftID) = data.KPI(shiftID,id);
        P{driverID}(shiftID,:) = data.KPI(shiftID, 7:end);
    end
end

%% NOW DOING PLOTS
markers = {'x', 'o', '*', '.', '+', 'diamond', 'square', '^', '<', '>', 'hexagram', 'pentagram', 'none', 'x'};
variables = ["$t_{pass}$", "$v_x$", "$a_x$", "$\omega$", "$\kappa_{road}$", "$d\kappa$"];

f = figure(1);
f.Position = [100 100 1800 550];
set(f,'defaulttextInterpreter','latex') ;
set(f, 'defaultAxesTickLabelInterpreter','latex');  
set(f, 'defaultLegendInterpreter','latex');

for i=1:size(P,2)
    subplot(2,5,i);
    surf([1:1:6], [1:1:10], (P{i}));
    ylabel("$np_i$"); xlabel("$P$");
    xticks([1:1:6]);yticks([1:1:10]);
    title(strcat("$P^{dr", num2str(i),"}$"));
    set(gca,'TickLabelInterpreter','latex');
    set(gca,'FontSize', 14);
end

%% FITTING ACCURACY

f = figure(2);
f.Position = [100 100 550 350];
set(f,'defaulttextInterpreter','latex') ;
set(f, 'defaultAxesTickLabelInterpreter','latex');  
set(f, 'defaultLegendInterpreter','latex');

surf([1:1:10], [1:1:10], NRMS');
ylabel("$np_i$"); xlabel("$dr_i$");
xticks([1:1:10]); yticks([1:1:10]);
if(id==2)
    title(strcat("$NRMS_{e_{\delta}}^{val}$"));
else
    title(strcat("$NRMS_{e_{\delta}}^{tr}$"));
end

set(gca,'TickLabelInterpreter','latex');
set(gca,'FontSize', 14);
