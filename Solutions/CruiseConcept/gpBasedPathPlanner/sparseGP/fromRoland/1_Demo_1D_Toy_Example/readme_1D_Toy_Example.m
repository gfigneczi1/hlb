%% This folder contains the code for the comparison of:
%    1)   SOGP (dictionary-based online GP, or the so-called sparse online
%         GP, widely used in most of the online GP papers): 
%         [1] Csató, L. and Opper, M., 2002. Sparse on-line Gaussian processes.  Neural computation, 14(3), pp.641-668. ---> the first paper!
%         [2] Chowdhary, G., Kingravi, H.A., How, J.P. and Vela, P.A., 2013, December. A Bayesian nonparametric approach to adaptive control using Gaussian processes. In 52nd IEEE Conference on Decision and Control (pp. 874-879). IEEE.
%    2)   ROGP (Recursive Online Sparse GP, proposed in Sec. 3) 
%    3)   DGP  (Dual GP, proposed n Sec. 4)

% and generates the figures in Sec. 5.1.


%% For SOGP case: 
% First, run main_dictionary_origin.m with "if_2nd_round = 0"
% The, run it again with "if_2nd_round = 1".

%% For ROGP case:
% First, run Case0_varGP.m (the result is save as "varGP.mat" in workspace)
% Second, run Case0_varGP_2nd.m (the GP update is based on the result of 1st iteration) 

%% For DGP case,:
% First, run Case1_DGP.m (the result is save as "res_DGP_1.mat" in workspace)
% second, run  Case1_DGP_2nd.m (the GP update is based on the result of 1st iteration) 

%% Routine for Dual GP
% 1. Data collection: Construct training data set D_0.
% 2. Offline training: Train the long-term GP with D_0, and initialize the short-term GP.
% 3. Online phase: Use the Dual GP for prediction for the task.
% 4. Offline retrain: Now we have new data from Step 3. Construct an augmented data set D_1, and go back to Step 2.
% 5. Online phase: Use the Dual GP for prediction for the repetitive/similar tasks. The experience of the past tasks are stored in the long-term GP.

%% Tunable parameters of Dual GP
% There are several parameters need to be tuned manually to get good predicition performance
% 1. lambda (forgetting factor)
% 2. Kernel hyperparameters of short-term GP
% 3. Initial posterior variance S0 for short-term GP
% 4. Initial inducing point: make it around the mean of the training input data to get good convergence of optimization


%% Trouble shooting
% In the ideal case, the estimation will conincide with its true value in the second iteration, as
% showin in good_result_DGP_2ndIter.fig.

% However, because of the optimzation function used in this code (Vargp toolbox
% written by Michalis Titsias), the convergence sometimes runs into a local
% minimum, and the regression results looks not that good. I think the
% reason is that there are several randomized initial parameters.

% In this case, you only need to run the Dual GP routine again. (Case1_DGP
% ---> Case1_DGP_2nd). Then you will see the expected result.


%% For comparison of the computational load:
% See the folder /MC_simulation
