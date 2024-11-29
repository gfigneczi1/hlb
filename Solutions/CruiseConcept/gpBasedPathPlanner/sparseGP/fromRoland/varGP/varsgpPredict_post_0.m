function [mustar, varstar] = varsgpPredict_post(model, Xtest)
%
%
persistent invKm L 

if isempty(invKm)

sigma2 = exp(2*model.Likelihood.logtheta);
if strcmp(model.indType, 'pseudoIns')
      Kmm = kernel(model.GP, model.Xu); 
      Knm = kernel(model.GP, model.X, model.Xu); 
      Kmn = Knm';
      Kmstar = kernel(model.GP, model.Xu, Xtest);
elseif strcmp(model.indType, 'weights')
      Kmm = kernelWeights(model);
      Knm = kernelWeights(model, model.X);
      Kmn = Knm';
      Kmstar = kernelWeights(model, Xtest)';
end

Kmm = Kmm+model.jitter*eye(model.m);
L   =  chol(Kmm)';  % K = LL';
invL = inv(L);
invKm = invL'*invL;

end
% V = L\Kmn;
% 
% D = sigma2;
% V = V./repmat(sqrt(D)',model.m,1);
% y = model.y./sqrt(D);
% 
% M = eye(model.m) + V*V'; % M = I+VD-1V'
% Lm = chol(M)';  % M = LmLm'

% make predictions 
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




