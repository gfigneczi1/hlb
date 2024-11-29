%% This is the code for the paper: Learning for Predictive Control: A Dual Gaussian Process Approach
% Before you start, add all the subfolders into the matlab path
% There are three demos prepared:
% 0. ROGP (Recursive Online Sparse GP) for 1D example
% 1. Comparison of three sparse online algorithms (SOGP, ROGP, DGP), corresponding the results in Sec. 5.1
% 2. Applications in predictive control: Comparison of Baseline MPC, LGPMPC, OGPMPC, and DGPMPC
%    Baseline MPC:  Naive nonlinear MPC without model augmentation
%    LGPMPC:        MPC strategy with offline trained GP
%    OGPMPC:        MPC strategy with offline trained GP which is recursively update online
%    DGPMPC:        MPC strategy with Dual GP structure

%% To start with the code
% Go to /0_Demo_ROGP_1D_Example, run Demo_Vargp_online.m. You will see how
% the ROGP updates its posterior under the condition that the target
% function is varying.

%% Next, go to  /1_Demo_1D_Toy_Example
%  Check the readme_1D_Toy_Example.m for more details.
%  Please note that this is a 1D  example. If you are using this code for
%  some real-world applications, please construct GP models for each dimension respectively, and train them seperately.
%  For example: Car-like scenario:              model_x & model_y
%               Quadcopter trajectory tracking: model_x & model_y & model_z

%% Third, go to /2-Demo_DGPMPC_Quadcopter
% This is a quadcopter trajectory tracking task, with time-varying
% external wind. The GP is offline trained under the condition of contant
% wind, then tested online with sin-like time-varying wind.

% Step 1. You can first run Case_1_baselineMPC.m see the naive result of baseline
% MPC without GP compensation.

% Step 2. Then, run Step_1_dataColle.m to collect the training data set. The
% quadcopter is forced to track a pseudo random trajectory with MPC.

% Step 3. Run Step_2_5_DGP_offlineTrain_SVGP.m to train the long-term GP
% offine and initialize the short-term GP.

% Step 4. After running Step_2_5_DGP_offlineTrain_SVGP.m, you can run
% Case_2_lGPMPC.m, or Case_3_OGPMPC.m, or Case_4_DGPMPC.m. 

% Step 5. After running Case_4_DGPMPC.m, now we have the "experience" one
% task. Then run Step_3_DGP_after_iteration_1st.m to train the long-term GP
% offline, and re-initialize the short term GP. 

% Step 6. Run Case_5_DGPMPC_2nd.m to the tracking performance in the second
% round.


% Please note that in this scenario we need to construct a seperate GP for
% each dimension. So the "mex" file is used to accelerate the online
% prediction.



%% Some important files for implementing the algorithm
%    ROSV_GP_lambda.m        ----> ROGP 
%    varsgpPredict_MM_x_c.m  ----> variational sparse GP with moment matching
%    DGP_varsgpPredict_x_c.m ----> dual GP with moment matching