function [mustar, varstar] = varsgpPredict_post_1(model, Xtest)
%
%

sigma2 = exp(2*model.Likelihood.logtheta);


Kmm = model.Kmm;
Kmn = model.Kmn;
Knm = model.Knm;
invKm = model.invKm;
L     = model.L;
Kmstar = kernel(model.GP, model.Xu, Xtest);


L   =  chol(Kmm)';  % K = LL';


% make predictions 
lst = L\Kmstar;

post_m = model.postm;
post_S = model.postS;

Kstar = kernel(model.GP, Xtest, [], 1); 

mustar   = Kmstar'*invKm*post_m;
varstar  = Kstar - sum(lst.^2,1)' + diag(Kmstar'*invKm*post_S*invKm*Kmstar);




