function evaluator_identifyControlModel(segments,config)

% This evaluator is responsible for revealing options how the control
% theoretic model can be used for modelling the steering behavior of the
% driver.

% Segmentor: segmentor_driverModel.

% @copyright (C) 2023 Robert Bosch GmbH.
% The reproduction, distribution and utilization of this file as
% well as the communication of its contents to others without express
% authorization is prohibited. Offenders will be held liable for the
% payment of damages. All rights reserved in the event of the grant
% of a patent, utility model or design.
% @version 1.0
% @date 2023-04-14

temp_folder_path = config.root;
plots_folder_name = 'plots';

for fileID=1:size(segments.segments,2)
    segment = segments.segments(fileID).segment;
    
    try
        
        save(fullfile(temp_folder_path, plots_folder_name, strcat(segments.segments(fileID).name,'.mat')),'segment');

        clear segment;
    catch e
        disp("The following file could not be evaluated:");
        disp(segments.segments(fileID).name);
        disp(e);
    end
end

end


