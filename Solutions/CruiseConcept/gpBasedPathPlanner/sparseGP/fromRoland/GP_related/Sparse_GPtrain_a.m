function [pu,qu,ru,hyp_1,hyp_2,hyp_3]  = Sparse_GPtrain_a(training_in_a,training_out_a,Nu)



N = length(training_in_a);

% randomly choose Nu inducing points as the initial input
index  = randsample(N, Nu);

attitude_u = training_in_a(index,:);  % this is the initial inducign points, Nu*D

[Nu,D] = size(attitude_u);


input    = training_in_a;

target_1 = training_out_a(:,1);
target_2 = training_out_a(:,2);
target_3 = training_out_a(:,3);

% for the first dimension

disp('GP training for p')

hyp_1 = [log(ones(D+1,1));log(1)];
para_1 = [reshape(attitude_u,Nu*D,1);hyp_1];

[para_1,~] = minimize(para_1, 'spgp_lik', -100, target_1,input,Nu);

pu = reshape(para_1(1:Nu*D,1),Nu,D);
hyp_1 = para_1(Nu*D+1:end,1);



% for the second dimension
disp('GP training for q')

hyp_2 = [log(ones(D+1,1));log(1)];
para_2 = [reshape(attitude_u,Nu*D,1);hyp_2];

[para_2,~] = minimize(para_2, 'spgp_lik', -100, target_2,input,Nu);

qu = reshape(para_2(1:Nu*D,1),Nu,D);
hyp_2 = para_2(Nu*D+1:end,1);



% for the third dimension
disp('GP training for r')

hyp_3 = [log(ones(D+1,1));log(1)];
para_3 = [reshape(attitude_u,Nu*D,1);hyp_3];

[para_3,~] = minimize(para_3, 'spgp_lik', -100, target_3,input,Nu);

ru = reshape(para_3(1:Nu*D,1),Nu,D);
hyp_3 = para_3(Nu*D+1:end,1);



% [mu1, var1] = gp(hyp_1, @infExact, [], covfunc, likfunc, input, target_1,input);
% [mu2, var2] = gp(hyp_2, @infExact, [], covfunc, likfunc, input, target_2,input);
% [mu3, var3] = gp(hyp_3, @infExact, [], covfunc, likfunc, input, target_3,input);
% 
% mu  =  [mu1 mu2 mu3];
% var =  [var1 var2 var3]; 
