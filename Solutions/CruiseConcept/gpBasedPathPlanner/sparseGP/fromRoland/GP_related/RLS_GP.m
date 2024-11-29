function [alpha_new,P_new,absErr,L] = RLS_GP(alpha,err,phi,P_old)

%% Input: alpha: weight 
%         err  : the error between the regression and true value
%         phi  : covariance matrix
P = P_old;

   lambda  = 0.96;
%   lambda  = 0.8;
% lambda  =  0.1;
% lambda  =  1;
% lambda  =  0.99;

tmp1    = lambda + ( phi * P * phi');

invTmp1 = 1/tmp1; 

L       = invTmp1 * P * phi';   % 20*1

tmpMat1 = P * (phi'*phi) * P;

P_new       = lambda^(-1)*(P -  invTmp1*tmpMat1);

alpha_new =  alpha + L*err;

absErr= abs(err);


end