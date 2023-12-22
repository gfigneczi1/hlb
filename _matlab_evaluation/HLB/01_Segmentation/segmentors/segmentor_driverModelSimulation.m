function [segment] = segmentor_driverModelSimulation(input_folder)
% Segmentation for Driver Model Simulation
% 1. Read in the input and output time series data for the Driver Model
% 2. Create a segment struct to save the following signals from simulation: 
% all inputs of the model (segment.input); trajectoryGlobal, 
% trajectoryPlannerFrame, corridorPlannerFrame, egoPoseGlobalFrame (segment.output)

    % Read in
    matFiles = dir(fullfile(input_folder,'/*.mat'));
    MatFilesTable = struct2table(matFiles);
    clear matFiles;
    sortedMatFilesTable = sortrows(MatFilesTable, 'date');
    matFiles = table2struct(sortedMatFilesTable);

    for fileID = 1:length(matFiles)
        if matFiles(fileID).name == "model_in.mat"
        	input = load(fullfile(matFiles(fileID).folder,matFiles(fileID).name));
        elseif matFiles(fileID).name == "model_out.mat"
            output = load(fullfile(matFiles(fileID).folder,matFiles(fileID).name));
        end
    end
    
    % Make sure signals are synchron
%     if input.model_in.GPS_time.Data == output.model_out.GPS_time.Data   
%     else
%         disp("Time synchronization problem with the inputs and outputs of the model")
%         return
%     end
    
    
    % Convert Input data
    segment.input.objectVelocityFront = (input.model_in.objectVelocityFront.Data);
    segment.input.objectDistanceFront = (input.model_in.objectDistanceFront.Data);
    segment.input.objectAccelerationFront = (input.model_in.objectAccelerationFront.Data);
    segment.input.GPS_time = (input.model_in.GPS_time.Data);
    segment.input.AccelerationX_ESP = (input.model_in.AccelerationX_ESP.Data);
    segment.input.AccelerationY_ESP = (input.model_in.AccelerationY_ESP.Data);
    segment.input.VelocityX_ESP = (input.model_in.VelocityX_ESP.Data);
    segment.input.yawRateESP = (input.model_in.yawRateESP.Data);
    segment.input.SteeringTorque = (input.model_in.SteeringTorque.Data);
    segment.input.SteeringAngle = (input.model_in.SteeringAngle.Data);
    segment.input.c01_right = (input.model_in.c01_right.Data);
    segment.input.c01_left = (input.model_in.c01_left.Data);
    segment.input.q_T0 = (input.model_in.objectVelocityFront.time);
    segment.input.c1 = (input.model_in.c1.Data);
    segment.input.c2 = (input.model_in.c2.Data);
    segment.input.X_abs = input.model_in.X_abs.Data;
    segment.input.Y_abs = input.model_in.Y_abs.Data;
    segment.input.theta_calc = input.model_in.theta_calc.Data;
    
    segment.input.start_x = 1432927.251124854 * ones(length(input.model_in.GPS_time.data),1);
    segment.input.start_y = 5232154.292243999 * ones(length(input.model_in.GPS_time.data),1);
    segment.input.GPS_status = 8 * ones(length(input.model_in.GPS_time.data),1);
    segment.input.KF_status = 1 * ones(length(input.model_in.GPS_time.data),1);
    segment.input.LaneChange_Approved = 0 * ones(length(input.model_in.GPS_time.data),1);
    segment.input.Right_Index = 0 * ones(length(input.model_in.GPS_time.data),1);
    segment.input.Left_Index = 0 * ones(length(input.model_in.GPS_time.data),1);
    segment.input.c02_left = 0 * ones(length(input.model_in.GPS_time.data),1);
    segment.input.c02_right = 0 * ones(length(input.model_in.GPS_time.data),1);

    
    % Convert Output data
    segment.output.trajectoryGlobalFrame = (output.model_out.trajectoryGlobalFrame.Data);
    
    % The following Outputs are commented out in the model as well and can
    % be added if needed
%     segment.output.trajectoryPlannerFrame = (output.model_out.trajectoryPlannerFrame.Data);
%     segment.output.corridorPlannerFrame.corridor_X = (output.model_out.corridorPlannerFrame.corridor_X.Data);
%     segment.output.corridorPlannerFrame.corridor_Y = (output.model_out.corridorPlannerFrame.corridor_Y.Data);
%     segment.output.corridorPlannerFrame.corridor_c1 = (output.model_out.corridorPlannerFrame.corridor_c1.Data);
%     segment.output.corridorPlannerFrame.corridor_c2 = (output.model_out.corridorPlannerFrame.corridor_c2.Data);
%     segment.output.corridorGlobalFrame.corridorGlobal_X = (output.model_out.corridorGlobalFrame.corridorGlobal_X.Data);
%     segment.output.corridorGlobalFrame.corridorGlobal_Y = (output.model_out.corridorGlobalFrame.corridorGlobal_Y.Data);
%     segment.output.corridorGlobalFrame.corridorGlobal_c1 = (output.model_out.corridorGlobalFrame.corridorGlobal_c1.Data);
%     segment.output.corridorGlobalFrame.corridorGlobal_c2 = (output.model_out.corridorGlobalFrame.corridorGlobal_c2.Data);
%     segment.output.egoPoseGlobalFrame = (output.model_out.egoPoseGlobalFrame.Data);
end