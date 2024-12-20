clear; close all; clc;

dataPathGp = "C:\git\KDP\publications\GP\results\fullGP";
dataPathSparseGp = "C:\git\KDP\publications\GP\results\sparseGP_fixedAcceleration\0_chosenParams";
dataPathPHTPM = "C:\git\KDP\publications\GP\results\PHTPM";
dataPathNARX = "C:\git\KDP\publications\GP\results\NARMAX\kpis"; % and ARX in the same files
dataPathLRM = "C:\git\KDP\publications\GP\results\LRM\kpis";
dataPathELDM = "C:\git\KDP\publications\GP\results\ELDM\kpi";

% loading the performance parameters
[NRMS_gp, NRMS_gp_val, RMS_gp, RMS_gp_val, trun_gp] = loadPerformanceFiles(dataPathGp,-1);
[NRMS_sgp, NRMS_sgp_val, RMS_sgp, RMS_sgp_val, trun_sgp] = loadPerformanceFiles(dataPathSparseGp,-1);
[NRMS_lrm, NRMS_lrm_val, RMS_lrm, RMS_lrm_val, trun_lrm] = loadPerformanceFiles(dataPathLRM,0);
[NRMS_eldm, NRMS_eldm_val, RMS_eldm, RMS_eldm_val, trun_eldm] = loadPerformanceFiles(dataPathELDM, 2);
[NRMS_ldm, NRMS_ldm_val, RMS_ldm, RMS_ldm_val, trun_ldm] = loadPerformanceFiles(dataPathELDM, 1);

[NRMS_arx, NRMS_arx_val, RMS_arx, RMS_arx_val, trun_arx] = loadPerformanceFiles(dataPathNARX, 3);
[NRMS_narx, NRMS_narx_val, RMS_narx, RMS_narx_val, trun_narx] = loadPerformanceFiles(dataPathNARX, 4);

[NRMS_phtpm, NRMS_phtpm_val, RMS_phtpm, RMS_phtpm_val, trun_phtpm] = loadPerformanceFiles(dataPathPHTPM, 5);

% merging to arrays
NRMS_sum = [mean(mean(NRMS_arx)) mean(mean(NRMS_narx)) mean(mean(NRMS_gp)) mean(mean(NRMS_sgp))  mean(mean(NRMS_lrm)) mean(mean(NRMS_ldm)) mean(mean(NRMS_eldm)) mean(mean(NRMS_phtpm)); ...
    median(median(NRMS_arx)) median(median(NRMS_narx)) median(median(NRMS_gp)) median(median(NRMS_sgp))  median(median(NRMS_lrm)) median(median(NRMS_ldm)) median(median(NRMS_eldm)) median(median(NRMS_phtpm)); ...
    (std(NRMS_arx)) (std(NRMS_narx)) (std(NRMS_gp)) std(std(NRMS_sgp))  std(std(NRMS_lrm)) std(std(NRMS_ldm)) std(std(NRMS_eldm)) (std(NRMS_phtpm)); ...
    mean(mean(NRMS_arx_val)) mean(mean(NRMS_narx_val)) mean(mean(NRMS_gp_val)) mean(mean(NRMS_sgp_val))  mean(mean(NRMS_lrm_val)) mean(mean(NRMS_ldm_val)) mean(mean(NRMS_eldm_val)) mean(mean(NRMS_phtpm_val)); ...
    median(median(NRMS_arx_val)) median(median(NRMS_narx_val)) median(median(NRMS_gp_val)) median(median(NRMS_sgp_val))  median(median(NRMS_lrm_val)) median(median(NRMS_ldm_val)) median(median(NRMS_eldm_val)) median(median(NRMS_phtpm_val)); ...
    (std(NRMS_arx_val)) (std(NRMS_narx_val)) (std(NRMS_gp_val)) std(std(NRMS_sgp_val))  std(std(NRMS_lrm_val)) std(std(NRMS_ldm_val)) std(std(NRMS_eldm_val)) (std(NRMS_phtpm_val)); ...
    (mean(trun_arx)) mean(mean(trun_narx)) mean(mean(trun_gp)) mean(mean(trun_sgp))  mean(mean(trun_lrm)) mean(mean(trun_ldm)) mean(mean(trun_eldm)) mean(mean(trun_phtpm))];

NRMS_sum(end,:) = NRMS_sum(end,:)/max(NRMS_sum(end,:));

f = spider_plot(NRMS_sum', "AxesLabels", {"$\overline{NRMS}_{e_{\delta}}^{train}$", "$NRMS_{e_{\delta}}^{train,med}$", "$\sigma_{e_{\delta}}^{train}$", "$\overline{NRMS}_{e_{\delta}}^{val}$", "$NRMS_{e_{\delta}}^{val,med}$", "$\sigma_{e_{\delta}}^{val}$", "$t_{run}$"});
set(f,'defaulttextInterpreter','latex') ;
set(f, 'defaultAxesTickLabelInterpreter','latex');  
set(f, 'defaultLegendInterpreter','latex');
set(gca,'TickLabelInterpreter','latex');

fprintf("NRMS mean values - train:\n\tNRMS_gp=%f\n\tNRMS_sgp=%f\n\tNRMS_lrm=%f\n\tNRMS_eldm=%f\n\tNRMS_ldm=%f\n\tNRMS_arx=%f\n\tNRMS_narx=%f\n\tNRMS_phtpm=%f", mean(mean(NRMS_gp)), ...
    mean(mean(NRMS_sgp)), mean(mean(NRMS_lrm)), mean(mean(NRMS_eldm)), mean(mean(NRMS_ldm)), mean(mean(NRMS_arx)), mean(mean(NRMS_narx)), mean(mean(NRMS_phtpm)));

fprintf("\nNRMS median values - train:\n\tNRMS_gp=%f\n\tNRMS_sgp=%f\n\tNRMS_lrm=%f\n\tNRMS_eldm=%f\n\tNRMS_ldm=%f\n\tNRMS_arx=%f\n\tNRMS_narx=%f\n\tNRMS_phtpm=%f", median(median(NRMS_gp)), ...
    median(median(NRMS_sgp)), median(median(NRMS_lrm)), median(median(NRMS_eldm)), median(median(NRMS_ldm)), median(median(NRMS_arx)), median(median(NRMS_narx)), median(median(NRMS_phtpm)));


fprintf("\n\nNRMS std values - train:\n\tNRMS_gp=%f\n\tNRMS_sgp=%f\n\tNRMS_lrm=%f\n\tNRMS_eldm=%f\n\tNRMS_ldm=%f\n\tNRMS_arx=%f\n\tNRMS_narx=%f\n\tNRMS_phtpm=%f", (std(NRMS_gp)), ...
    std(std(NRMS_sgp)), std(std(NRMS_lrm)), std(std(NRMS_eldm)), std(std(NRMS_ldm)), (std(NRMS_arx)), (std(NRMS_narx)), (std(NRMS_phtpm)));

fprintf("\n\nNRMS mean values - validation:\n\tNRMS_gp=%f\n\tNRMS_sgp=%f\n\tNRMS_lrm=%f\n\tNRMS_eldm=%f\n\tNRMS_ldm=%f\n\tNRMS_arx=%f\n\tNRMS_narx=%f\n\tNRMS_phtpm=%f", mean(mean(NRMS_gp_val)), ...
    mean(mean(NRMS_sgp_val)), mean(mean(NRMS_lrm_val)), mean(mean(NRMS_eldm_val)), mean(mean(NRMS_ldm_val)), mean(mean(NRMS_arx_val)), mean(mean(NRMS_narx_val)), mean(mean(NRMS_phtpm_val)));

fprintf("\n\nNRMS std values - validation:\n\tNRMS_gp=%f\n\tNRMS_sgp=%f\n\tNRMS_lrm=%f\n\tNRMS_eldm=%f\n\tNRMS_ldm=%f\n\tNRMS_arx=%f\n\tNRMS_narx=%f\n\tNRMS_phtpm=%f\n", (std(NRMS_gp_val)), ...
    std(std(NRMS_sgp_val)), std(std(NRMS_lrm_val)), std(std(NRMS_eldm_val)), std(std(NRMS_ldm_val)), (std(NRMS_arx_val)), (std(NRMS_narx_val)), (std(NRMS_phtpm_val)));

function [NRMS, NRMS_val, RMS, RMS_val, trun] = loadPerformanceFiles(data, modifier)
    KPIFiles = dir(fullfile(data,"KPI*.mat"));
    for i=1:length(KPIFiles)
        a = load(fullfile(KPIFiles(i).folder, KPIFiles(i).name));
        if (iscell(a.KPI))
            numberOfShifts = size(a.KPI,2);
        else
            numberOfShifts = 1;
        end
        for shiftId = 1:numberOfShifts
            if (modifier==-1)
                NRMS(i,shiftId) = a.KPI{shiftId}(end-2);
                NRMS_val(i,shiftId) = a.KPI{shiftId}(2);
                RMS(i,shiftId) = a.KPI{shiftId}(end-4);
                RMS_val(i,shiftId) = a.KPI{shiftId}(1);
                trun(i,shiftId) = a.KPI{shiftId}(end);
            elseif (modifier == 0)
                NRMS(i,shiftId) = a.KPI{shiftId}(6);
                NRMS_val(i,shiftId) = a.KPI{shiftId}(2);
                RMS(i,shiftId) = a.KPI{shiftId}(5);
                RMS_val(i,shiftId) = a.KPI{shiftId}(1);
                trun(i,shiftId) = a.KPI{shiftId}(end);
            elseif (modifier == 1)
                NRMS(i,shiftId) = a.KPI{shiftId}(6);
                NRMS_val(i,shiftId) = a.KPI{shiftId}(2);
                RMS(i,shiftId) = a.KPI{shiftId}(5);
                RMS_val(i,shiftId) = a.KPI{shiftId}(1);
                trun(i,shiftId) = a.KPI{shiftId}(size(a.KPI{shiftId},2)/2);
            elseif (modifier == 2)
                NRMS(i,shiftId) = a.KPI{shiftId}(size(a.KPI{shiftId},2)/2+6);
                NRMS_val(i,shiftId) = a.KPI{shiftId}(size(a.KPI{shiftId},2)/2+2);
                RMS(i,shiftId) = a.KPI{shiftId}(size(a.KPI{shiftId},2)/2+5);
                RMS_val(i,shiftId) = a.KPI{shiftId}(size(a.KPI{shiftId},2)/2+1);
                trun(i,shiftId) = a.KPI{shiftId}(end);
            elseif (modifier == 3)
                NRMS(i,shiftId) = a.KPI(1,2);
                NRMS_val(i,shiftId) = a.KPI(1,5);
                RMS(i,shiftId) = a.KPI(1,1);
                RMS_val(i,shiftId) = a.KPI(1,4);
                trun(i,shiftId) = a.KPI(1,end);
            elseif (modifier == 4)
                NRMS(i,shiftId) = a.KPI(2,1);
                NRMS_val(i,shiftId) = a.KPI(2,5);
                RMS(i,shiftId) = a.KPI(2,1);
                RMS_val(i,shiftId) = a.KPI(2,4);
                trun(i,shiftId) = a.KPI(2,end);
            elseif (modifier == 5)
                NRMS(i,shiftId) = a.KPI(1,5);
                NRMS_val(i,shiftId) = a.KPI(1,2);
                RMS(i,shiftId) = a.KPI(1,4);
                RMS_val(i,shiftId) = a.KPI(1,1);
                trun(i,shiftId) = a.KPI(1,end);
            end
        end
    end

    % removing driver 2
    NRMS = [NRMS(1:5,:); NRMS(7:end,:)];
    NRMS_val = [NRMS(1:5,:); NRMS(7:end,:)];
    RMS = [NRMS(1:5,:); NRMS(7:end,:)];
    RMS_val = [NRMS(1:5,:); NRMS(7:end,:)];
    trun = [trun(1:5,:); trun(7:end,:)];
end