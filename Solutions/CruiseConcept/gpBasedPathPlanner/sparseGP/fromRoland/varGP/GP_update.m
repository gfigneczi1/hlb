%% Update the GP weight alpha in a recursive manner
% based on the offline trained Sparse GP

function [Alpha,err] = GP_update(state_hist, statedot_hist, R_3_hist, F_hist,torque_hist,Alpha,mu_p,mu_a,hyp_1,hyp_2,hyp_3,hyp_4,hyp_5,hyp_6)

global BV_1 BV_2 BV_3 BV_4 BV_5 BV_6 P10 P20 P30 P40 P50 P60

alpha_p1 = Alpha(:,1);
alpha_p2 = Alpha(:,2);
alpha_p3 = Alpha(:,3);
alpha_a1 = Alpha(:,4);
alpha_a2 = Alpha(:,5);
alpha_a3 = Alpha(:,6);

%% For position

target_p = getDiffOnline_p(state_hist(end,:),statedot_hist(end,:),R_3_hist(end,:),F_hist(end));   % 1*3

err_p = target_p - mu_p;

% err_p = mu_p(2)-kern([state_hist(end,4:9),F_hist(end)],BV_2,hyp_2)*alpha;

phi_1 = kern([state_hist(end,4:9),F_hist(end)],BV_1,hyp_1);  % 1*20
phi_2 = kern([state_hist(end,4:9),F_hist(end)],BV_2,hyp_2);
phi_3 = kern([state_hist(end,4:9),F_hist(end)],BV_3,hyp_3);

[alpha_p1,P1,err_p1] = RLS_GP(alpha_p1,err_p(1),phi_1,P10);
[alpha_p2,P2,err_p2] = RLS_GP(alpha_p2,err_p(2),phi_2,P20);
[alpha_p3,P3,err_p3] = RLS_GP(alpha_p3,err_p(3),phi_3,P30);

P10 = P1;
P20 = P2;
P30 = P3;


%% For attitude

target_a = getDiffOnline_a(state_hist(end,:),statedot_hist(end,:),torque_hist(end,:));   % 1*3
err_a = target_a - mu_a;

phi_4 = kern([state_hist(end,10:12),torque_hist(end,:)],BV_4,hyp_4);
phi_5 = kern([state_hist(end,10:12),torque_hist(end,:)],BV_5,hyp_5);
phi_6 = kern([state_hist(end,10:12),torque_hist(end,:)],BV_6,hyp_6);

[alpha_a1,P4,err_a1] = RLS_GP(alpha_a1,err_a(1),phi_4,P40);
[alpha_a2,P5,err_a2] = RLS_GP(alpha_a2,err_a(2),phi_5,P50);
[alpha_a3,P6,err_a3] = RLS_GP(alpha_a3,err_a(3),phi_6,P60);
P40 = P4;
P50 = P5;
P60 = P6;

Alpha = [alpha_p1 alpha_p2 alpha_p3 alpha_a1 alpha_a2 alpha_a3];
% err = [err_p1 err_p2 err_p3 err_a1 err_a2 err_a3];
err  = [err_p err_a];
end

