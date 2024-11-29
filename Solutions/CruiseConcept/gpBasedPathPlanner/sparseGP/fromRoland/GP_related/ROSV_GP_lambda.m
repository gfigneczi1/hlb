function model_1 = ROSV_GP_lambda(model,err,phi,lambda)

%% Input: model 
%         err  : the error between the regression and true value
%         phi  : covariance matrix

    model_1 = model;
  
%     lambda = 0.985;
%     lambda = 0.98;
%     lambda = 1;

    S = model.postS;
    m = model.postm;
%     S0 = model.postS0;

%     lambda = trace(S)/trace(S0);

    Gk = lambda + phi*S*phi';

    iGk = 1/Gk;

    Lk   = iGk * S * phi';   %M*1

    Hk = S * phi' * iGk;

    S_new = lambda^(-1)*(S -  Hk*Gk*Hk');


    model_1.postS = S_new;
    model_1.postm = m + Lk*err;




end