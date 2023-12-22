function matlab_evaluation_cover_test(root)
% ========================  
% #########################################################################
% Great thanks to Balint Wirnhardt, the original creator of this file.
% Although it's been almost totally changed... But anyway, thanks!

% The following layers are presented in the script:
% - Input reading
% - Segmentation
% - Evaluation

% This is the test version of matlab_evaluation_cover.m!

    clc;
    % Input reading
    if (nargin == 0)
        root = fullfile('..','..','_temp');
    end
           
    if exist(root,'dir')
        config = jsondecode(fileread(fullfile(root,'config.json')));
    end
    config.root = root;
    
    configuration_segmentors = {"MapValidation", "MapBand", "DriverModel", "DriverModel", "DriverModel", "DriverModel"};
    configuration_evaluators = {"MapValidation", "MapBand", "NodePointDefinition", "RoadGeneration", "DriverModelAnalysis", "DriverModelLearning"};
    
    for i=1:length(configuration_segmentors)
        clc;
        config.segmentation_profile = configuration_segmentors{i};
        config.evaluation_profile = configuration_evaluators{i};
        try
            % Segmentation
            [segments, config] = segmentor(config,root);
            % Evaluation
            evaluator(segments, config);
            resultString(i) = (strcat('OK:',{' '},config.segmentation_profile,{' '},config.evaluation_profile));
        catch e
            resultString(i) = (strcat('FAILED:',{' '},config.segmentation_profile,{' '},config.evaluation_profile, {' '}, 'with error:', e));
        end
    end
    disp(resultString);

end

