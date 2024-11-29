function [par_new,absErr] = ROGP_SVI(err,phi,par)

%% Input: alpha: weight 
%         err  : the error between the regression and true value
%         phi  : covariance matrix
lambda = 0.95;

Gk = lambda + phi*par.S*phi';

iGk = 1/Gk;

Lk   = iGk * par.S * phi';   % 20*1

Hk = par.S * phi' * iGk;

S_new = lambda^(-1)*(par.S -  Hk*Gk*Hk');


par.S = S_new;
par.m = par.m + Lk*err;

absErr= abs(err);
par_new = par;


end