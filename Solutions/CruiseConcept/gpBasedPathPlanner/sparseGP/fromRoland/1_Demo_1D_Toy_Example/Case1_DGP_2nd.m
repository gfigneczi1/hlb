%% run this after Case1_DGP.m
% do not clear the work space 

clear;
% rand('seed', 1e5);

load DGP_1st

flag_update =1;
Xtrain_new = xn;
Ytrain_new = ff;

Xtrain =  Xtrain_new';
Ytrain = Ytrain_new';

sigma = 0.2;
 Ytrain =  Ytrain + sigma*randn(length(Ytrain),1);  

[nn,D] = size(Xtrain);
%% batchtrainThelong-termGP
% type of likelihood 
options.Likelihood = 'Gaussian';

%number of inducing variables.
options.m = 20; 

% type of inducing variables
options.indType = 'pseudoIns'; % 'pseudoIns' or 'weights'
options.objectFunc = 'var';

new_model = varsgpCreate('seard', Xtrain, Ytrain, options);

Xuinit = mean(Xtrain(:)) + 2*randn(new_model.m,1); 
new_model.Xu = Xuinit;
new_model.Xuinit = Xuinit;

trops(1) = 300; % number of iterations
trops(2) = 1; % type of optimizations
trops(3) = 1; 

% initialization of the model hyperparameters 
logtheta0(1:D,1) = log((max(Xtrain)-min(Xtrain))'/2);
logtheta0(D+1,1) = 0.5*log(var(Ytrain,1)); 
logtheta0(D+2,1) = 0.5*log(var(Ytrain,1)/4);  
new_model.GP.logtheta = logtheta0(1:end-1);
new_model.Likelihood.logtheta = logtheta0(end);

%% offline train
% train the model by optimizing over the kernel hyperparameters  
% the inducing variables parameters 
[new_model margLogL] = varsgpTrain(new_model, trops);

N = length(Xtrain);
for i = 1:N

    [mu_x(i,:),var_x(i,:)] = varsgpPredict_post_x(new_model,Xtrain(i));
 

end
std_x  = sqrt(var_x);

%% re-initialize the short term GP

new_Smodel = varsgpCreate('seard', Xtrain,Ytrain, options);
index  = randsample(nn,new_Smodel.m);  
Xuinit =  Xtrain(index,:); 
new_Smodel.Xu = Xuinit;

logtheta0(1:D,1) = 2*ones(D,1); % smoother estimation
logtheta0(D+1,1) = -2;
logtheta0(D+2,1) = -5;
new_Smodel.GP.logtheta = logtheta0(1:end-1);
new_Smodel.Likelihood.logtheta = logtheta0(end);

new_Smodel.postm = zeros(new_Smodel.m,1);
new_Smodel.postS = eye(new_Smodel.m)*1e0;
[~,~,Knm,invKm,L,Kmm] = computeKminv(new_Smodel,Xtrain);
new_Smodel.Kmm = Kmm;
new_Smodel.Kmn = Knm';
new_Smodel.Knm = Knm;
new_Smodel.L   = L;
new_Smodel.invKm = invKm;


Xtest = xn';

%% simulation
  s1 = tic;

for i  = 1:N
    
%    [mu,s2] = varsgpPredict_post_1(model, Xtest(i));
   [mu,s2] = DGP_varsgpPredict_x(new_model,new_Smodel,Xtest(i),0);
   
    err(i) = ff(i)-mu;
    
     if  flag_update == 1

           lambda = 0.999;
      [phix,~,~,~,~] = computeKminv(new_Smodel, Xtest(i));  % new kernel
         
      Smodel = ROSV_GP_lambda(new_Smodel, err(i),phix,lambda);
         
%
    end
   



mustar(i,:) = mu;
varstar(i,:) = s2;

   
   
end
mustar_DGP_2 = mustar;
varstar_DGP_2 = varstar;
save('res_DGP_2.mat','mustar_DGP_2','varstar_DGP_2');


toc(s1)

varstar(1:3,:) = zeros(3,1);
figure;
hold on;
fillColor = [236 241 249]/255;
hold on; grid on;box on;
plot(Xtrain,Ytrain,'k.','linewidth',1);
set(gca,'Fontsize',15);
plot(Xtest, mustar, 'b-.', 'LineWidth', 2) % mean predictions in blue
plot(Xtest, ff', 'r', 'LineWidth', 2) % 
set(gca, 'fontsize', 18);
box on;grid on;
legend('95% Conf','Training data','Mean','True');