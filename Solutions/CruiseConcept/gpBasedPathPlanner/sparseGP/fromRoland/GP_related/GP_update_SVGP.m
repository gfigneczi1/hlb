function [model_x1,model_y1,model_z1,err] = GP_update_SVGP(test_in,target,mu,model_x,model_y,model_z,lambda)



err = target - mu;

[phi_1,~,~,~,~] = computeKminv(model_x, test_in);  % new kernel  1xM
[phi_2,~,~,~,~] = computeKminv(model_y, test_in);  % new kernel
[phi_3,~,~,~,~] = computeKminv(model_z, test_in);  % new kernel


model_x1 = ROSV_GP_lambda(model_x,err(1),phi_1,lambda);
model_y1 = ROSV_GP_lambda(model_y,err(2),phi_2,lambda);
model_z1 = ROSV_GP_lambda(model_z,err(3),phi_3,lambda);

end
