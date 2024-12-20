function params = loadParameters(modelID, DriverID)

    %% HERE GIVE THE PATH TO THE PARAMETER FILES IN THE FUNCTION OF THE MODEL ID
    switch modelID
        case "gp"
            pathToParams = "C:\git\KDP\publications\GP\results\fullGP";
        case "sgp"
            pathToParams = "C:\git\KDP\publications\GP\results\sparseGP_fixedAcceleration\0_chosenParams";
        case "phtpm"
            pathToParams = "C:\git\KDP\publications\GP\results\PHTPM";
        case "ldm"
            pathToParams = "C:\git\KDP\publications\GP\results\ELDM\kpi";
        case "eldm"
            pathToParams = "C:\git\KDP\publications\GP\results\ELDM\kpi";
        case "arx"
            pathToParams = "C:\git\KDP\publications\GP\results\NARMAX\kpis";
        case "narx"
            pathToParams = "C:\git\KDP\publications\GP\results\NARMAX\kpis";
        case "lrm"
            pathToParams = "C:\git\KDP\publications\GP\results\LRM\kpis";
    end
    
    paramFiles = dir(fullfile(pathToParams,"ETA_*"));

    %% HERE CALL THE RIGHT PARAMETER READ IN FILE WHICH WILL RETURN SPECIFIED STRUCT OF THE PARAMETERS
    % This consists of two different parts:
    % 1: common fields: c_in/s_in/c_out/s_out (the necessary data for
    % normalization and centralization, and other specific fields, such as
    % hyp_opt_array (for GP), net (for DNN based solutions), P (for
    % LDM/ELDM/LRM)...etc

    for fileID = 1:length(paramFiles)
        paramDriverID = str2num(paramFiles(fileID).name(strfind(paramFiles(fileID).name, 'driver_')+7:strfind(paramFiles(fileID).name, '.')-1));
        if (paramDriverID == DriverID)
            paramData = load(fullfile(paramFiles(fileID).folder, paramFiles(fileID).name));

            % The specific fields
            switch modelID
                case "gp"
                    for npID=1:length(paramData.ETA)
                        params.hyp_opt_array{npID} = paramData.ETA(npID).hyp_opt;
                        params.output_estimation{npID} = paramData.ETA(npID).output_estimation;
                        params.input_estimation{npID} = paramData.ETA(npID).input_estimation;
                    end
                case "sgp"
                    for npID=1:length(paramData.ETA)
                        params.hyp_opt_array{npID} = paramData.ETA(npID).hyp_opt;
                        params.output_estimation{npID} = paramData.ETA(npID).output_estimation;
                        params.input_estimation{npID} = paramData.ETA(npID).input_estimation;
                    end
                case "phtpm"
                    params.net = paramData.ETA.net;
                    params.XTrain = paramData.ETA.XTrain;
                    params.TTrain = paramData.ETA.TTrain;
                    params.XVal = paramData.ETA.XVal;
                    params.TVal = paramData.ETA.TVal;
                case "ldm"
                    params.P = paramData.ETA.P;
                case "eldm"
                    params.P_left = paramData.ETA.P_left;
                    params.P_right = paramData.ETA.P_right;
                case "lrm"
                    params.P = paramData.ETA.P;
                case "narx"
                    params.sysNARX = paramData.ETA.sysNARX;
                case "arx"
                    params.sysARX = paramData.ETA.sysARX;
            end

            % The common fields
            for npID = 1:length(paramData.ETA)
                params.c_in = paramData.ETA(npID).normFactors.c_in; % common for all GPs
                params.s_in = paramData.ETA(npID).normFactors.s_in; % common for all GPs

                % workaround for LDM-ELDM sharing the input data
                if (modelID == "eldm")
                    params.c_in = [params.c_in params.c_in];
                    params.s_in = [params.s_in params.s_in];
                end
                
                params.c_out = paramData.ETA(npID).normFactors.c_out; % one at each node point
                params.s_out = paramData.ETA(npID).normFactors.s_out; % one at each node point

                % workaround for PHTPM containing the output as input, too
                if (modelID == "phtpm")
                    params.c_in = [params.c_in params.c_out];
                    params.s_in = [params.s_in params.s_out];
                    params.s_out(1:15) = params.s_out;
                    params.c_out(1:15) = params.c_out;
                end
            end
            break;
        else
            params = [];
        end
    end

end

