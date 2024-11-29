function [mu_MM,var_MM,var_long, var_short, var_cross] = DGP_varsgpPredict_x(model,Smodel,m_in,S_in)
% for single dimension



[n,D]   =   size(model.X);


sigma2f= exp(2*model.GP.logtheta(D+1));
sigma2f_s = exp(2*Smodel.GP.logtheta(D+1));

%% long-term part


Kmm = model.Kmm;
Kmn = model.Kmn;
Knm = model.Knm;

invKm = model.invKm;



%% short-term part


Kmm_s = Smodel.Kmm;
Kmn_s = Smodel.Kmn;
Knm_s = Smodel.Knm;

invKm_s = Smodel.invKm;


post_m = model.postm;
post_S = model.postS;
Spost_m = Smodel.postm;
Spost_S = Smodel.postS;



B = invKm -invKm*post_S*invKm;
B_s = invKm_s -invKm_s*Spost_S*invKm_s;


%% MM 

% long-term mean
inp  = bsxfun(@minus,model.Xu,m_in);

beta = invKm*post_m;
iLambda = diag(exp(-2*model.GP.logtheta(1:D)));     % DxD, inversed squared length scale
R    = S_in + diag(exp(2*model.GP.logtheta(1:D)));
iR   = iLambda * (eye(D) - (eye(D)+S_in*iLambda)\(S_in*iLambda));
T    = inp*iR;
c    = sigma2f/sqrt(det(R))*exp(sum(model.GP.logtheta(1:D)));
q    = c*exp(-sum(T.*inp,2)/2);
qb   = q.*beta;

% short-term mean
inp_s  = bsxfun(@minus,Smodel.Xu,m_in);

beta_s = invKm_s*Spost_m;
iLambda_s = diag(exp(-2*Smodel.GP.logtheta(1:D)));     % DxD, inversed squared length scale
R_s    = S_in + diag(exp(2*Smodel.GP.logtheta(1:D)));
iR_s   = iLambda_s * (eye(D) - (eye(D)+S_in*iLambda_s)\(S_in*iLambda_s));
T_s    = inp_s*iR_s;
c_s    = sigma2f_s/sqrt(det(R_s))*exp(sum(Smodel.GP.logtheta(1:D)));
q_s    = c_s*exp(-sum(T_s.*inp_s,2)/2);
qb_s   = q_s.*beta_s;

% DGP predict mean
mu_MM = sum(qb) +  sum(qb_s) ;


% long-term  var
v    = bsxfun(@rdivide,inp,exp(model.GP.logtheta(1:D))');   % xu-mu/sqrt(iLambda)
log_k = 2*model.GP.logtheta(D+1) - sum(v.*v,2)/2;

Zeta_i  = bsxfun(@rdivide,inp,exp(2*model.GP.logtheta(1:D)'));   % DxD

R = S_in*diag(exp(-2*model.GP.logtheta(1:D))+exp(-2*model.GP.logtheta(1:D))) + eye(D);
t = 1./sqrt(det(R));

Q = t*exp(bsxfun(@plus,log_k,log_k') + maha(Zeta_i,-Zeta_i, R\S_in/2));

A = beta*beta';
A = A - B;
A = A.*Q;

var_long = sum(sum(A));
var_long = var_long + sigma2f;

% short-term var
v_s    = bsxfun(@rdivide,inp_s,exp(Smodel.GP.logtheta(1:D))');   % xu-mu/sqrt(iLambda)
log_ks = 2*Smodel.GP.logtheta(D+1) - sum(v_s.*v_s,2)/2;

Zeta_is  = bsxfun(@rdivide,inp_s,exp(2*Smodel.GP.logtheta(1:D)'));   % DxD

R_s = S_in*diag(exp(-2*Smodel.GP.logtheta(1:D))+exp(-2*Smodel.GP.logtheta(1:D))) + eye(D);
t_s = 1./sqrt(det(R_s));

Q_s = t_s*exp(bsxfun(@plus,log_ks,log_ks') + maha(Zeta_is,-Zeta_is, R_s\S_in/2));

A_s = beta_s*beta_s';
A_s= A_s - B_s;
A_s = A_s.*Q_s;

var_short = sum(sum(A_s));
var_short = var_short + sigma2f_s;

% cross var
R_x = S_in*diag(exp(-2*model.GP.logtheta(1:D))+exp(-2*Smodel.GP.logtheta(1:D))) + eye(D);
t_x = 1./sqrt(det(R_x));

Q_x = t_x*exp(bsxfun(@plus,log_k,log_ks') + maha(Zeta_i,-Zeta_is, R_x\S_in/2));

A_x = beta*beta_s';
A_x = A_x.*Q_x;
var_cross = 2*sum(sum(A_x));

% total var = long + short + cross
var_MM = var_long + var_short + var_cross - mu_MM*mu_MM';





