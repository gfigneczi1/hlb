close all;
clear;

matFiles = dir(fullfile("C:\git\KDP\publications\GP\results\kernel_verification_setparams\KPI_input*.mat"));
id = 2; %6: train data, 2: validation data
selectedKernel = "LIN-PP3";

for i=1:length(matFiles)
    data = load(fullfile(matFiles(i).folder, matFiles(i).name));
    driverID = str2num(matFiles(i).name(strfind(matFiles(i).name, 'driver_')+7:strfind(matFiles(i).name, '.')-1));
    inputID = str2num(matFiles(i).name(strfind(matFiles(i).name, 'input_')+6:strfind(matFiles(i).name, '_driver_')-1));
    for shiftID = 1:numel(data.KPI)
        NRMS(driverID, inputID, shiftID) = data.KPI{shiftID}(id);
        eta{driverID, inputID}(shiftID,:) = data.KPI{shiftID}(9:end);
        sigma(driverID, inputID, shiftID) = exp(data.KPI{shiftID}(end));
    end
end

%% NOW DOING PLOTS FOR ALL KERNELS
for driverID=1:size(NRMS,1)
    for inputID=1:size(NRMS,2)
        NRMS_avgd(driverID,inputID) = mean(NRMS(driverID,inputID,:));
    end
end

markers = {'x', 'o', '*', '.', '+', 'diamond', 'square', '^', '<', '>', 'hexagram', 'pentagram', 'none', 'x'};
variables = ["$k_1$", "$k_2$", "$k_3$", "$k_4$", "$k_5$"];

f = figure(3);
f.Position = [100 100 450 200];
set(f,'defaulttextInterpreter','latex') ;
set(f, 'defaultAxesTickLabelInterpreter','latex');  
set(f, 'defaultLegendInterpreter','latex');

boxplot(NRMS_avgd(:,:));
xticklabels(variables);
grid on;
ylim([0,0.12]);

set(gca,'FontSize', 14);
if (id==2)
    title("Estimation accuracy - validation data");
    ylabel("$NRMS_{e_\delta}^{val}$");
else
    title("Estimation accuracy - train data");
    ylabel("$NRMS_{e_\delta}^{tr}$");
end

set(gca,'TickLabelInterpreter','latex');

switch selectedKernel
    case "k_8"
        lambda_SE = [1:1:6];
        sigma_SE = 7;
        lamdba_PP = 8:1:13;
        sigma_PP = 14;
        k = 8;
    case "LIN-PP3"
        lambda_LIN = [1:1:8];
        lambda_PP = 9:1:16;
        sigma_PP = 17;
        k = 4;
    otherwise
end

if (selectedKernel =="SE-PP3")
    % SE-ARD
    f = figure(1);
    f.Position = [100 100 1800 550];
    set(f,'defaulttextInterpreter','latex') ;
    set(f, 'defaultAxesTickLabelInterpreter','latex');  
    set(f, 'defaultLegendInterpreter','latex');
    
    for i=1:size(eta,1)
        subplot(2,4,i);
        surf([1:1:6], [1:1:10], exp(eta{i,k}(:,lambda_SE)));
        ylabel("$GP_i$"); xlabel("$\lambda_{SE}$");
        xticks([1:1:6]);yticks([1:1:10]);
        title(strcat("$\lambda_{SE}^{dr", num2str(i),"}$"));
        set(gca,'TickLabelInterpreter','latex');
        set(gca,'FontSize', 14);
    end
    
    % PP3-ARD
    f = figure(2);
    f.Position = [100 100 1800 550];
    set(f,'defaulttextInterpreter','latex') ;
    set(f, 'defaultAxesTickLabelInterpreter','latex');  
    set(f, 'defaultLegendInterpreter','latex');
    
    for i=1:size(eta,1)
        subplot(2,4,i);
        surf([1:1:6], [1:1:10], exp(eta{i,k}(:,lamdba_PP)));
        ylabel("$GP_i$"); xlabel("$\lambda_{P3}$");
        xticks([1:1:12]);yticks([1:1:10]);
        %zlim([-0,3]);
        title(strcat("$\lambda_{P3}^{dr", num2str(i),"}$"));
        set(gca,'TickLabelInterpreter','latex');
        set(gca,'FontSize', 14);
    end
elseif (selectedKernel =="LIN-PP3")
    % LIN-ARD
    f = figure(1);
    f.Position = [100 100 1800 550];
    set(f,'defaulttextInterpreter','latex') ;
    set(f, 'defaultAxesTickLabelInterpreter','latex');  
    set(f, 'defaultLegendInterpreter','latex');
    
    for i=1:size(eta,1)
        subplot(2,4,i);
        surf([1:1:8], [1:1:10], exp(eta{i,k}(:,lambda_LIN)));
        ylabel("$GP_i$"); xlabel("$\lambda_{LIN}$");
        xticks([1:1:8]);yticks([1:1:10]);
        title(strcat("$\lambda_{LIN}^{dr", num2str(i),"}$"));
        set(gca,'TickLabelInterpreter','latex');
        set(gca,'FontSize', 14);
    end
    
    % PP3-ARD
    f = figure(2);
    f.Position = [100 100 1800 550];
    set(f,'defaulttextInterpreter','latex') ;
    set(f, 'defaultAxesTickLabelInterpreter','latex');  
    set(f, 'defaultLegendInterpreter','latex');
    
    for i=1:size(eta,1)
        subplot(2,4,i);
        surf([1:1:8], [1:1:10], exp(eta{i,k}(:,lambda_PP)));
        ylabel("$GP_i$"); xlabel("$\lambda_{P3}$");
        xticks([1:1:8]);yticks([1:1:10]);
        %zlim([-0,3]);
        title(strcat("$\lambda_{P3}^{dr", num2str(i),"}$"));
        set(gca,'TickLabelInterpreter','latex');
        set(gca,'FontSize', 14);
    end
    
    % PP3 sigma
    for i=1:size(eta,1)
        sigmas_PP3(:,i) = exp(eta{i,k}(:,sigma_PP));
    end
end