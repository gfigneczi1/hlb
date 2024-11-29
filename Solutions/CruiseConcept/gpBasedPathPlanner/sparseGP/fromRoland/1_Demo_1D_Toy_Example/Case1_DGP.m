clc;clear; close all

addpath('misc');
addpath('SPGP_dist');
addpath('docs');

flag_update =1;
dT = 0.05;

minX = -5;
maxX = 0;           % training set
minX_Test = -5;
maxX_Test = 10;      % testing set
% 


fx = @(x) 2*sin(x)+0.05*x.^2+2;
fg = @(x) fx(x)-2*cos(x)-0.2*x;



xn = minX_Test:dT:maxX_Test;  % column
N = length(xn);
xn1 = xn(1:ceil(N/2)); xn2 = xn(ceil(N/2)+1:end);
f1 = fx(xn1);
f0 = fx(xn);
f2 = fg(xn2);
% f2 = fx(xn2);
ff = [f1 f2];


sigma = 0.2;
% choose the measurement data among xn as the training set

Nm = 200;    % training set size  （big）
xm = sort(minX + rand(1,Nm)*(maxX-minX));  % choose randmoly
ym = fx(xm)+ sigma*randn(1,Nm);  % noisy data 
Xtrain = xm';
Ytrain = ym';
xn = minX_Test:dT:maxX_Test;  % column'
Xtest = xn';

if flag_update == 1
   Ytest = fg(xn)';
else 
    Ytest = fx(xn)';
end
    
 


my = 0;
Ytrain = Ytrain - my;
[n D] = size(Xtrain);

% type of likelihood 
options.Likelihood = 'Gaussian';

%number of inducing variables.
options.m = 20; 

% type of inducing variables
options.indType = 'pseudoIns'; % 'pseudoIns' or 'weights'
options.objectFunc = 'var';

model = varsgpCreate('seard', Xtrain, Ytrain, options);

index  = randsample(Nm,model.m);  
Xuinit =  Xtrain(index,:); 
model.Xu = Xuinit;

trops(1) = 100; % number of iterations
trops(2) = 1; % type of optimizations
trops(3) = 1; 

% initialization of the model hyperparameters 
logtheta0(1:D,1) = log((max(Xtrain)-min(Xtrain))'/2);
logtheta0(D+1,1) = 0.5*log(var(Ytrain,1)); 
logtheta0(D+2,1) = 0.5*log(var(Ytrain,1)/4);  
model.GP.logtheta = logtheta0(1:end-1);
model.Likelihood.logtheta = logtheta0(end);

%% offline train for long-term-GP

[model margLogL] = varsgpTrain(model, trops);

% initialize the short-term GP 
Smodel = varsgpCreate('seard', Xtrain,Ytrain, options);
index  = randsample(Nm,Smodel.m);  
Xuinit =  Xtrain(index,:); 
Smodel.Xu = Xuinit;

logtheta0(1:D,1) = 2*ones(D,1); % smoother estimation
logtheta0(D+1,1) = -1;
logtheta0(D+2,1) = -5;
Smodel.GP.logtheta = logtheta0(1:end-1);
Smodel.Likelihood.logtheta = logtheta0(end);

Smodel.postm = zeros(Smodel.m,1);
Smodel.postS = eye(Smodel.m)*1e0;
[~,~,Knm,invKm,L,Kmm] = computeKminv(Smodel,Xtrain);
Smodel.Kmm = Kmm;
Smodel.Kmn = Knm';
Smodel.Knm = Knm;
Smodel.L   = L;
Smodel.invKm = invKm;

 lambda = 0.97;
%% Simulation starts

  s1 = tic;

for i  = 1:N
    

   [mu,s2,var_long, var_short, var_cross] = DGP_varsgpPredict_x(model,Smodel,Xtest(i),0); % make prediction: long/short 
   
    err(i) = ff(i)-mu;
    
     if  flag_update == 1

          
      [phix,~,~,~,~] = computeKminv(Smodel, Xtest(i));  % new kernel
         
      Smodel = ROSV_GP_lambda(Smodel, err(i),phix,lambda);       % at the same time, update the short-term GP recursively
         
%
    end
  

    mustar(i,:) = mu;
    varstar(i,:) = s2;

   
   
end

  toc(s1)


mustar_DGP_1 = mustar;
varstar_DGP_1 = varstar;
save('res_DGP_1.mat','mustar_DGP_1','varstar_DGP_1');

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


save('DGP_1st.mat','xn','ff');