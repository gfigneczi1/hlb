function [mustar, varstar] = varsgpPredict_post_y(model, Xtest)
%
%
sigma2 = exp(2*model.Likelihood.logtheta);


Kmm = model.Kmm;
Kmn = model.Kmn;
Knm = model.Knm;
Kmstar = kernel(model.GP, model.Xu, Xtest);

Kmm = Kmm + model.jitter*eye(model.m);
L   =  chol(Kmm)';  % K = LL';
invL = inv(L);
invKm = invL'*invL;





% make predictions 


Kmstar = kernel(model.GP, model.Xu, Xtest);
lst = L\Kmstar;
% MinvL = Lm*Lm'/L; 
% LmL = Lm\L';

% post_m = MinvL\(V*y);
% post_S = LmL'*LmL;
 post_m = model.postm;
 post_S = model.postS;

Kstar = kernel(model.GP, Xtest, [], 1); 

mustar   = Kmstar'*invKm*post_m;

varstar  = Kstar - sum(lst.^2,1)' + diag(Kmstar'*invKm*post_S*invKm*Kmstar);







