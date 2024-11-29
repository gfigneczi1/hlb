function [A,Kmstar,Knm,invKm,L,Kmm] = computeKminv(model, Xtest)
%COMPUTEKNMKMMINV [A,Knm,Kmminv,Lmm,Kmm] = computeKnmKmminv(covfunc,loghyp,x,z)
%   Compute the term A = k(x,z)k(z,z)^{-1} which is commonly encountered.
%
Kmm = kernel(model.GP, model.Xu)+model.jitter*eye(model.m);
Knm = kernel(model.GP, model.X, model.Xu); 
L   =  chol(Kmm)';  % K = LL';
invL = inv(L);
invKm = invL'*invL;

Kmstar = kernel(model.GP, model.Xu, Xtest); % Km*

A = Kmstar'*invKm;

% Kmm = feval(covfunc, loghyp, z) + 1e-10*eye(size(z,1));
% Lmm = jit_chol(Kmm);
% Kmminv = invChol(Lmm);
% Knm = feval(covfunc, loghyp, x, z);
% A = Knm*Kmminv;
end