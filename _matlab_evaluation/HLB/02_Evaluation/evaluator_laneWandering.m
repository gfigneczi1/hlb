function evaluator_laneWandering(segments,config)

% This evaluator is responsible for generating various plots to analyse
% lane wandering performance of different drivers.
% The following KPI-s are proposed to be used:
% - straight line:
% -- compensation points (where the steering torque is applied to compensate stray effect)
% -- distribution of lane offset and CLP
% -- oncoming traffic? follow traffic?

% Segmentor: segmentor_driverModel.

% @copyright (C) 2022 Robert Bosch GmbH.
% The reproduction, distribution and utilization of this file as
% well as the communication of its contents to others without express
% authorization is prohibited. Offenders will be held liable for the
% payment of damages. All rights reserved in the event of the grant
% of a patent, utility model or design.
% @version 1.0
% @date 2023-09-04


%% Step 1: 
% Snippeting with different frequencies to define the right cut off
% frequency. This function does not return anything, but saves plots and
% data mat files into temp/plots/frequencies subfolder. Then,
% standaloneScripts/laneWandering_offsetErrorSeparation must be run to
% collect the files of given drivers and calculate the cost of each
% frequencies.
%functional_laneWanderingFrequencyDefinition(segments, config);

%% Step 2:
% intervention point calculation. This function does not return anything
% but calculates the distribution of the offset error extrema assoicated
% with the intervention point and saves the corresponding plot.
functional_laneWanderingInterventionPoints(segments, config);

%% Step 3: 
% optimizitation of the mpc params on the snippets OR re-simulate 
% optimized parameters.
functional_laneWanderingMpc(segments, config);

end

