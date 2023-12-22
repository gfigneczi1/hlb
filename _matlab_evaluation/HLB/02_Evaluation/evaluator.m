function evaluator(segments, config)
% MAIN evlauator function which calls the subsequent evaluators based on
% the configuration.

% @copyright (C) 2022 Robert Bosch GmbH.
% The reproduction, distribution and utilization of this file as
% well as the communication of its contents to others without express
% authorization is prohibited. Offenders will be held liable for the
% payment of damages. All rights reserved in the event of the grant
% of a patent, utility model or design.
% @version 1.0
disp(strcat("Evaluation",{' '},config.evaluation_profile,{' '},"started"));

% just for testing
temp_folder_path = config.root;
plots_folder_name = 'plots';
    
%create or reinitialize _temp\plots folder
if ~isfolder(temp_folder_path)
    mkdir(temp_folder_path);
end
if ~isfolder(fullfile(temp_folder_path, plots_folder_name))
    mkdir(fullfile(temp_folder_path, plots_folder_name));
else
    rmdir(fullfile(temp_folder_path, plots_folder_name),'s');
    mkdir(fullfile(temp_folder_path, plots_folder_name));
end
switch config.evaluation_profile
    case "MapValidation"
        evaluator_mapValidation(segments,config);
    case "MapBand"
        evaluator_mapBands(segments, config);
    case "NodePointDefinition"
        evaluator_nodePointDefinition(segments,config);
    case "RoadGeneration"
        evaluator_roadGeneration(segments,config);
    case "DriverModelAnalysis"
        evaluator_driverModelAnalysis(segments,config);
    case "FullDriverModel"
        evaluator_fullDriverModel(segments,config);
    case "DriverModelLearning"
        evaluator_driverModelLearning(segments,config);
    case "MapValidationJkk"
        evaluator_mapValidation_jkk(segments,config);
    case "CompareModel"
        evaluator_compareDriverModel(segments,config);
    case "LDM"
        evaluator_olLDMResim(segments, config);
    case "ControlModel"
        evaluator_identifyControlModel (segments, config);
    case "LaneWandering"
        evaluator_laneWandering(segments,config);
    case "OptimalPath"
        evaluator_optimalPath(segments,config);
        
end