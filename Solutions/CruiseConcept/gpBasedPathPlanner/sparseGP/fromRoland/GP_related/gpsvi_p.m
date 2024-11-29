function [elbo_1 elbo_2 elbo_3,par_1,par_2,par_3]  = gpsvi_p(training_in_p,training_out_p,test_in,Nu,cf1,cf2,cf3)
%% offline train: get the q(u|y)=N(u|m,s)


input    = training_in_p;

target_1 = training_out_p(:,1);
target_2 = training_out_p(:,2);
target_3 = training_out_p(:,3);

%% for each dimension

disp('GP training for X')
[elbo_1,par_1] = svi_learn(input,target_1,test_in,Nu,cf1,[]);  % column input

disp('GP training for Y')
[elbo_2,par_2] = svi_learn(input,target_2,test_in,Nu,cf2,[]); 

disp('GP training for Z')
[elbo_3,par_3] = svi_learn(input,target_3,test_in,Nu,cf3,[]); 

%  elbo = [elbo_1 elbo_2 elbo_3];

end