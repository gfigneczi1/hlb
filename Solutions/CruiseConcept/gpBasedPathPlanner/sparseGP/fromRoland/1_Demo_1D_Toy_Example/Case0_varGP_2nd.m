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
% fx = @(x) 2*sin(x)+0.05*x.^2+2;
% % fx = @(x) 2*sin(x)+0.05*x+2
% fg = @(x) fx(x)-0.5*cos(x)-0.1*x.^2-0.2*x;

fx = @(x) 2*sin(x)+0.05*x.^2+2;
% fg = @(x) fx(x)-0.5*cos(x)-0.2*x;
fg = @(x) fx(x)-2*cos(x)-0.2*x;



xn = minX_Test:dT:maxX_Test;  % column
N = length(xn);
xn1 = xn(1:ceil(N/2)); xn2 = xn(ceil(N/2)+1:end);
f1 = fx(xn1);
f0 = fx(xn);
f2 = fg(xn2);
% f2 = fx(xn2);
ff = [f1 f2];
% fx = @(x) 2*sin(x)+0.05*x.^2+2;
% xn = minX:0.01:maxX;  % column

% figure;
% plot(xn,ff); hold on 
% plot(xn,f0); legend('f1','f0')

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





load varGP
% model.postS = eye(model.m)*1e2;   % !!!

lambda  = 0.99;  % forgetting factor
%% online prediction 

s1 = tic;

for i  = 1:N
    
   [mu,s2] = varsgpPredict_post_1(model, Xtest(i));
   
    err(i) = ff(i)-mu;
    
     if  flag_update == 1
         t0 = tic;
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
 
        %% new
%         Gk   =  Kstar - phi*Kmstar + phi*S*phi' + sigma;
%         iGk  = inv(Gk);    
%         Lk   = S*phi'*iGk;
%         S_new = (S - Lk*phi*S);
% 
        model.postS = S_new;
        model.postm = m + Lk*err(i);
        lambda_hist(i) = lambda;

        s0(i) = toc(t0);
    end
   
   
mustar(i,:) = mu;
varstar(i,:) = s2;

   
   
end

t1= toc(s1)
mustar(1,:) = mustar(2,:); % avoid some numerical issues
varstar(1:3,:) = zeros(3,1);

mustar_varGP_2 = mustar;
varstar_varGP_2 = varstar;
save('res_varGP_2.mat','mustar_varGP_2','varstar_varGP_2');
%% if postS is positive definite or not
% try chol(PostS_hist(:,:,1000))
%     disp('Matrix is symmetric positive definite.')
% catch ME
%     disp('Matrix is not symmetric positive definite')
% end
% 
% d = eig(PostS_hist(:,:,1000));
% isposdef = all(d > 0)


% add noise variance sigma2 to get a prediction for the ys

% varstar = varstar + exp(2*model.Likelihood.logtheta);

% display prediction in test data
% % % run firstly the full GP to know the ground truth 
% covfunc = {'covSum', {'covSEard','covNoise'}};
% covfunc = @covSEard;
% likfunc = @likGauss;
% l=1; sf=1; sigma = 0.2;
% 
% hyp     = struct('cov',log([l sf]),'lik',log(sigma));
% hyp     = minimize(hyp, 'gp', -100, @infGaussLik, [], covfunc, likfunc, Xtrain, Ytrain);
% 
% [muFull s2Full] = gp(hyp, @infGaussLik, [], covfunc, likfunc, Xtrain, Ytrain, Xtest);

% [logthetaFullGP, fullF] = minimize(logtheta0, 'gp', -200, covfunc, Xtrain, Ytrain);
% [muFull s2Full] = gpr(logthetaFullGP, covfunc, Xtrain, Ytrain, Xtest);


% varstar(1:5,:) = zeros(5,1);

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
% plot(Xuinit,-3*ones(size(Xuinit)),'r+','linewidth',2,'markersize',10);
% axis([min(Xtest) max(Xtest) -3 3]);
set(gca, 'fontsize', 18);
%set(gca, 'YTick', []);
%set(gca, 'XTick', []);
box on;grid on;
legend('95% Conf','Training data','Mean','True');

% std = sqrt(s2Full);
% figure('Name','Full GP');
% patch([Xtest' fliplr(Xtest)'],[muFull-2*std;flipud(muFull+2*std)]',[236 241 249]/255, 'EdgeColor', 'none');
% hold on; grid on;box on;
% plot(Xtest,muFull,'linewidth',1.5);
% plot(Xtrain,Ytrain,'r.','linewidth',1.5);
% plot(Xtest,f(Xtest),'linewidth',1.5);
% legend('95% Conf','Mean','Training data','True')
% title('Full GP')
% set(gca,'Fontsize',15);

% figure;
% plot(mustar,'b-','linewidth',2);hold on
% plot(ff','r:','linewidth',2);

% figure;
% plot(varstar);
% legend('Predicted variance');

% figure;
% subplot(3,1,1);
% plot(Xtest,m_hist,'b-','linewidth',2);
% legend('norm m');
% set(gca,'fontsize',16,'Fontname','times new roman');
% 
% subplot(3,1,2);
% plot(Xtest,S_hist,'b-','linewidth',2);
% legend('norm S');
% set(gca,'fontsize',16,'Fontname','times new roman');
% 
% subplot(3,1,3);
% plot(Xtest,lambda_hist,'b-','linewidth',2);
% legend('\lambda');
% set(gca,'fontsize',16,'Fontname','times new roman');