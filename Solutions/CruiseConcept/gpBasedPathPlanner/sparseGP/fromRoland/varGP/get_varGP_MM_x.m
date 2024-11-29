clc;clear;

load gpmodel
model = model_x;

syms m_in [1 7]
syms S_in [7 7]
% syms lambda

x = [m_in1 m_in2 m_in3 m_in4 m_in5 m_in6]';
u = m_in(7);

[n,D]   =   size(model.X);
sigma2 = exp(2*model.Likelihood.logtheta);
sigma2f = exp(2*model.GP.logtheta(D+1));

Kmm = model.Kmm;
Kmn = model.Kmn;
Knm = model.Knm;
invKm =  model.invKm;
L = model.L;
Xu = model.Xu;

 post_m = model.postm;
 post_S = model.postS;


 for i = 1: length(Xu)

    inp(i,:) = Xu(i,:)-m_in;
 end

 beta = invKm*post_m;
iLambda = diag(exp(-2*model.GP.logtheta(1:D)));     % DxD, inversed squared length scale
R    = S_in + diag(exp(2*model.GP.logtheta(1:D)));
iR   = iLambda * (eye(D) - (eye(D)+S_in*iLambda)\(S_in*iLambda));
T    = inp*iR;
c    = sigma2f/sqrt(det(R))*exp(sum(model.GP.logtheta(1:D)));
q    = c*exp(-sum(T.*inp,2)/2);
qb   = q.*beta;
mu_MM = sum(qb);

v    = bsxfun(@rdivide,inp,exp(model.GP.logtheta(1:D))');   % xu-mu/sqrt(iLambda)
log_k = 2*model.GP.logtheta(D+1) - sum(v.*v,2)/2;

Zeta_i  = bsxfun(@rdivide,inp,exp(2*model.GP.logtheta(1:D)'));   % DxD

R = S_in*diag(exp(-2*model.GP.logtheta(1:D))+exp(-2*model.GP.logtheta(1:D))) + eye(D);
t = 1./sqrt(det(R));

Q = t*exp(bsxfun(@plus,log_k,log_k') + maha(Zeta_i,-Zeta_i, R\S_in/2));

A = beta*beta';
A = A - B;
A = A.*Q;

var_MM = sum(sum(A));
var_MM = var_MM + sigma2f;

var_MM = var_MM - mu_MM*mu_MM';


dmu_dx =  matlabFunction(mu,'File','dmu1_dx','Vars',{m_in,S_in});






