function [segments, config] = segmentor(config, input_folder)
disp(strcat("Segmentation",{' '},config.segmentation_profile,{' '},"started"));
switch config.segmentation_profile
    case "MapValidation"
        segments = segmentor_mapValidation(input_folder);
    case "MapBand"
        segments = segmentor_map_band(input_folder);
    case "DriverModel"
        segments = segmentor_driverModel(input_folder);
    case "Simulation"
        segments = segmentor_driverModelSimulation(input_folder);
    case "MapValidationJkk"
        segments = segmentor_mapValidation_jkk(input_folder);
    case "Correction"
        % TO BE obsoleted!
        segments = segmentor_oldMeasCorrector(input_folder, config);
    case "LDM"
        segments = segmentor_olLDMResim(input_folder);
end