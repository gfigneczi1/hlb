%% using the posterior mean and variance to make predictions 
% origin from /working on


clc;clear; close all
%  randn('seed', 1e5);
% rand('seed', 1e5);
addpath('misc');
addpath('SPGP_dist');
addpath('docs');
% rng('shuffle');

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
    
 

% figure;
% plot(Xtrain,Ytrain,'.');

%   my = mean(Ytrain); 
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
logtheta0(1:D,1) =2*ones(D,1); % smoother estimation
logtheta0(D+1,1) = -1;
logtheta0(D+2,1) = -5;
model.GP.logtheta = logtheta0(1:end-1);
model.Likelihood.logtheta = logtheta0(end);

%% offline train


 model.postm = zeros(model.m,1);

 model.postS = eye(model.m)*1e2;

[~,~,Knm,invKm,L,Kmm] = computeKminv(model,Xtrain);
model.Kmm = Kmm;
model.Kmn = Knm';
model.Knm = Knm;
model.L   = L;
model.invKm = invKm;





lambda  = 0.97;  % forgetting factor
%% online prediction 

s1 = tic;

for i  = 1:N
    
   [mu,s2] = varsgpPredict_post_1(model, Xtest(i));
   
    err(i) = ff(i)-mu;
    
     if  flag_update == 1
         
         [phi,Kmstar,~,~,~] = computeKminv(model, Xtest(i));  % new kernel,  1xM
         
%          [model,Lk] = ROSV_GP_lambda(model,err(i),phi);
       

        S = model.postS;
        m = model.postm;

%     lambda = trace(S)/trace(S0);
    
%         % old
        Gk = lambda + phi*S*phi';

        iGk = 1/Gk;

        Lk  = iGk * S * phi';   % 20*1

        Hk = S * phi' * iGk;

%          Qd =  (Sd*phi'*phi*Sd)/(1+phi*Sd*phi');
       Qd = 0;

        S_new = lambda^(-1)*(S -  Hk*Gk*Hk'+ Qd) ;
 

% 
        model.postS = S_new;
        model.postm = m + Lk*err(i);
        lambda_hist(i) = lambda;
    end
   
   
mustar(i,:) = mu;
varstar(i,:) = s2;
m_hist(i,:) = norm(model.postm);
S_hist(i,:) = trace(model.postS);
PostS_hist(:,:,i)= model.postS;
   
   
end
mustar(1,:) = mustar(2,:); % avoid some numerical issues
varstar(1:3,:) = zeros(3,1);

t1= toc(s1)

mustar_varGP_1 = mustar;
varstar_varGP_1 = varstar;
save('varGP.mat','model');
save('res_varGP_1.mat','mustar_varGP_1','varstar_varGP_1','Xtest','ff');


figure;
hold on;
fillColor = [236 241 249]/255;
fill([Xtest; Xtest(end:-1:1)], [mustar;mustar(end:-1:1)]...
               + 2*[sqrt(varstar); -sqrt(varstar(end:-1:1))], fillColor,'EdgeColor','none');
hold on; grid on;box on;
% plot(Xtest,muFull,'linewidth',1.5);
%  plot(Xtrain,Ytrain,'k.','linewidth',1);
% plot(Xtest,fx(Xtest)-my,'linewidth',1.5);
set(gca,'Fontsize',15);

      
% plot(Xtrain, Ytrain,'.y','markersize',5); % data points in magenta

plot(Xtest, mustar, 'b-.', 'LineWidth', 2) % mean predictions in blue
plot(Xtest, ff', 'r', 'LineWidth', 2) % 
plot(Xtest, mustar + 2*sqrt(varstar),'-','LineWidth',1,'color',fillColor) % plus/minus 2 std deviation
plot(Xtest, mustar - 2*sqrt(varstar),'-','LineWidth',1,'color',fillColor) % plus/minus 2 std deviation

plot(model.Xu,-3*ones(size(model.Xu)),'b+','linewidth',2,'markersize',10);
plot(Xuinit,-3*ones(size(Xuinit)),'r+','linewidth',2,'markersize',10);
% axis([min(Xtest) max(Xtest) -3 3]);
set(gca, 'fontsize', 18);

box on;grid on;
legend('95% Conf','Training data','Mean','True');

