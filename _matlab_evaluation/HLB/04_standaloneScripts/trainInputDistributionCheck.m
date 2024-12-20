clear;close all; clc;

pathToParams = "C:\git\KDP\publications\GP\results\sparseGP_fixedAcceleration\0_chosenParams";
paramGP = dir(fullfile(pathToParams,"ETA_*"));

for driverID = 1:13
    for fileID = 1:length(paramGP)
        paramDriverID = str2num(paramGP(fileID).name(strfind(paramGP(fileID).name, 'driver_')+7:strfind(paramGP(fileID).name, '.')-1));
        if (paramDriverID == driverID)
            paramData = load(fullfile(paramGP(fileID).folder, paramGP(fileID).name));
            for npID=1:10
                GP_params.hyp_opt_array{npID} = paramData.ETA(npID).hyp_opt;
                GP_params.output_estimation{npID} = paramData.ETA(npID).output_estimation;
                GP_params.input_estimation{npID} = paramData.ETA(npID).input_estimation;
            end
            if (isfield(paramData.ETA(npID).normFactors, 'c_in'))
                GP_params.c_in = paramData.ETA(npID).normFactors.c_in; % common for all GPs
                GP_params.s_in = paramData.ETA(npID).normFactors.s_in; % common for all GPs
                GP_params.c_out(1:10) = paramData.ETA(npID).normFactors.c_out; % one at each node point
                GP_params.s_out = paramData.ETA(npID).normFactors.s_out; % one at each node point
            else
                GP_params.c_in = paramData.ETA(npID).normFactors(1:7); % common for all GPs
                GP_params.s_in = paramData.ETA(npID).normFactors(8:14); % common for all GPs
                GP_params.c_out(1:10) = paramData.ETA(npID).normFactors(15); % one at each node point
                GP_params.s_out = paramData.ETA(npID).normFactors(16:25); % one at each node point
            end
            break;
        else
            paramData = [];
        end
    end

    if (~isempty(paramData))
        % now loop through all inputs and generation 2D histogram from each
        m = size(GP_params.input_estimation{1}, 2);
        for i=1:size(GP_params.input_estimation{1}, 2)
            for j=1:size(GP_params.input_estimation{1}, 2)
                subplot(m, m, (i-1)*m+j);
                [counts, edges] = histcounts2(GP_params.input_estimation{1}(:,i), GP_params.input_estimation{1}(:,j), 51);
                imagesc(edges(1), edges(2), counts)
                colormap(jet)  % choose a colormap
                colorbar  % add a colorbar to the plot
                axis xy  % set the axis orientation to x-y
            end
        end
    end
end