%% Update the GP weight alpha in a recursive manner
% based on the offline trained Sparse GP

function [gp_new_trans,err] = GP_update(gpmodel_trans,mu_p,target_p,x,u)

global  P10 P20 P30 i

gp_new_trans = gpmodel_trans;
alpha_p1 = gpmodel_trans.alpha(:,1);
alpha_p2 = gpmodel_trans.alpha(:,2);
alpha_p3 = gpmodel_trans.alpha(:,3);


%% For position

err_p = target_p' - mu_p;

% err_p = mu_p(2)-kern([state_hist(end,4:9),F_hist(end)],BV_2,hyp_2)*alpha;
newstate = [x(7:9,i);x(4:6,i);u(3,i)]';  % 1*7

phi_1 = get_phi(gpmodel_trans,newstate,1);
phi_2 = get_phi(gpmodel_trans,newstate,2);
phi_3 = get_phi(gpmodel_trans,newstate,3);

[alpha_p1,P1,err_p1] = RLS_GP(alpha_p1,err_p(1),phi_1,P10);
[alpha_p2,P2,err_p2] = RLS_GP(alpha_p2,err_p(2),phi_2,P20);
[alpha_p3,P3,err_p3] = RLS_GP(alpha_p3,err_p(3),phi_3,P30);

P10 = P1;
P20 = P2;
P30 = P3;

err = [err_p1,err_p2,err_p3];

gp_new_trans.alpha(:,1) = alpha_p1;
gp_new_trans.alpha(:,2) = alpha_p2;
gp_new_trans.alpha(:,3) = alpha_p3;

end

