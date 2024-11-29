function [xu,yu,zu,hyp_1,hyp_2,hyp_3]  = Sparse_GPtrain_P(training_in_p,training_out_p,Nu)
% training_input_p: Nu*7      % are the pseudo inputs
% training_output_p: Nu*3     % corresponding targets


% likfunc = @likGauss;
% covfunc = @covSEard;

N = length(training_in_p);

% randomly choose Nu inducing points as the initial input
index  = randsample(N, Nu);

position_u = training_in_p(index,:);  % this is the initial inducign points, Nu*D

[Nu,D] = size(position_u);

input    = training_in_p;

target_1 = training_out_p(:,1);
target_2 = training_out_p(:,2);
target_3 = training_out_p(:,3);

%% for the first dimension

disp('GP training for X')
% hyp_1=struct('cov',ones(D+1,1),'lik',log(1),'sd',[],'ell',[],'K',[],'Ky',[]);

hyp_1 = [log(ones(D+1,1));log(1)];
para_1 = [reshape(position_u,Nu*D,1);hyp_1];

[para_1,~] = minimize(para_1, 'spgp_lik', -100, target_1,input,Nu);

xu = reshape(para_1(1:Nu*D,1),Nu,D);
hyp_1 = para_1(Nu*D+1:end,1);


%% for the second dimension
disp('GP training for Y')

hyp_2=[log(ones(D+1,1));log(1)];
para_2 = [reshape(position_u,Nu*D,1);hyp_2];

[para_2,~] = minimize(para_2, 'spgp_lik', -100, target_2,input,Nu);

yu = reshape(para_2(1:Nu*D,1),Nu,D);
hyp_2 = para_2(Nu*D+1:end,1);

%% for the third dimension
disp('GP training for Z')

hyp_3=[log(ones(D+1,1));log(1)];
para_3 = [reshape(position_u,Nu*D,1);hyp_3];

[para_3,~] = minimize(para_3, 'spgp_lik', -100, target_3,input,Nu);

zu = reshape(para_3(1:Nu*D,1),Nu,D);
hyp_3 = para_3(Nu*D+1:end,1);



