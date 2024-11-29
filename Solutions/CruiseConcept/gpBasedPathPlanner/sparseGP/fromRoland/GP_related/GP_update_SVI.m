%% Update the GP weight alpha in a recursive manner
% based on the offline trained Sparse GP

function [par_1_new,par_2_new,par_3_new,err] = GP_update_SVI(state_hist, statedot_hist, R_3_hist, F_hist,mu_p,cf1,cf2,cf3,par_1,par_2,par_3)



%% For position

target_p = getDiffOnline_p(state_hist(end,:),statedot_hist(end,:),R_3_hist(end,:),F_hist(end));   % 1*3

err_p = target_p - mu_p;

% err_p = mu_p(2)-kern([state_hist(end,4:9),F_hist(end)],BV_2,hyp_2)*alpha;

[phi_1,~,~,~,~] = computeKnmKmminv(cf1.covfunc,par_1.loghyp,[state_hist(end,4:9),F_hist(end)],par_1.z);  % new kernel
[phi_2,~,~,~,~] = computeKnmKmminv(cf2.covfunc,par_2.loghyp,[state_hist(end,4:9),F_hist(end)],par_2.z);  % new kernel
[phi_3,~,~,~,~] = computeKnmKmminv(cf3.covfunc,par_3.loghyp,[state_hist(end,4:9),F_hist(end)],par_3.z);  % new kernel
% 
% phi_1 = kern([state_hist(end,4:9),F_hist(end)],BV_1,hyp_1);  % 1*20
% phi_2 = kern([state_hist(end,4:9),F_hist(end)],BV_2,hyp_2);
% phi_3 = kern([state_hist(end,4:9),F_hist(end)],BV_3,hyp_3);

[par_1_new,err_p1] = ROGP_SVI(err_p(1),phi_1,par_1);
[par_2_new,err_p2] = ROGP_SVI(err_p(2),phi_2,par_2);
[par_3_new,err_p3] = ROGP_SVI(err_p(3),phi_3,par_3);


err = [err_p1,err_p2,err_p3];

end

