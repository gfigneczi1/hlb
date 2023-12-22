function matlab_evaluation_cover(root)
% ========================  
% #########################################################################
% Great thanks to Balint Wirnhardt, the original creator of this file.
% Although it's been almost totally changed... But anyway, thanks!

% The following layers are presented in the script:
% - Input reading
% - Segmentation
% - Evaluation

    % Input reading
    if (nargin == 0)
        root = fullfile('..','..','_temp');
    end
           
    if exist(root,'dir')
        config = jsondecode(fileread(fullfile(root,'config.json')));
    end
    config.root = root;
    
    % Segmentation
    [segments, config] = segmentor(config,root);
    
    % Evaluation
    evaluator(segments, config);
    disp('Evaluation succeeded');

end

