%% Code for Recursive online sparse GP 
% Training data: [-5,5] --> for offline training
% Test data: [5,10]     --> for online testing
% The target function varies from fx to fg at x =7.5



clc;clear; close all
randn('seed', 1e5);
rand('seed', 1e5);

flag_update = 1;   % 0 --> only sparse variational GP; 1--> recursive online sparse GP

minX = -5;
maxX = 5;           % training set
minX_Test = -5;
maxX_Test = 10;      % testing set

fx = @(x) 2*sin(x)+0.05*x.^2+2;
fg = @(x) fx(x)-0.5*cos(x)-0.1*x.^2-0.2*x;

xn = minX_Test:0.01:maxX_Test;  % column
N = length(xn);
xn1 = xn(1:ceil(N/2)); xn2 = xn(ceil(N/2)+1:end);
f1 = fx(xn1);
f0 = fx(xn);
f2 = fg(xn2);
ff = [f1 f2];

sigma = 0.2;
% choose the measurement data among xn as the training set

Nm = 500;    % training set size  （big）
xm = sort(minX + rand(1,Nm)*(maxX-minX));  % choose randmoly
ym = fx(xm)+ sigma*randn(1,Nm);  % noisy data 
Xtrain = xm';
Ytrain = ym';
xn = minX_Test:0.01:maxX_Test;  % column'
Xtest = xn';

if flag_update == 1
   Ytest = fg(xn)';
else 
    Ytest = fx(xn)';
end
    
 

%  my = mean(Ytrain); 
my = 0;
Ytrain = Ytrain - my;
[n D] = size(Xtrain);

% type of likelihood 
options.Likelihood = 'Gaussian';

%number of inducing variables.
options.m = 10; 

% type of inducing variables
options.indType = 'pseudoIns'; % 'pseudoIns' or 'weights'
options.objectFunc = 'var';

model = varsgpCreate('se', Xtrain, Ytrain, options);
% Fix seeds
% randn('seed', 2e5);
% rand('seed', 2e5);
Xuinit = mean(Xtrain(:)) + 0.5*randn(model.m,1); 
model.Xu = Xuinit;

trops(1) = 1000; % max number of iterations
trops(2) = 1; % type of optimizations
trops(3) = 1; 

% initialization of the model hyperparameters 
logtheta0(1:D,1) = log((max(Xtrain)-min(Xtrain))'/2);
logtheta0(D+1,1) = 0.5*log(var(Ytrain,1)); 
logtheta0(D+2,1) = 0.5*log(var(Ytrain,1)/4);  
model.GP.logtheta = logtheta0(1:end-1);
model.Likelihood.logtheta = logtheta0(end);

%% offline train
% train the model by optimizing over the kernel hyperparameters  
% the inducing variables parameters 
[model margLogL] = varsgpTrain(model, trops);
margLogL(end)

if flag_update == 1
    model.postS = eye(10)*1e0;
end

lambda  = 0.995;  % forgetting factor
%% online prediction 
for i  = 1:N
    
   [mu,s2] = varsgpPredict_post_1(model, Xtest(i));
   
    err(i) = ff(i)-mu;
    
     if  flag_update == 1
         
         [phi,~,~,~,~] = computeKminv(model, Xtest(i));  % new kernel
         
          model = ROSV_GP_lambda(model,err(i),phi,lambda);   % main function for the ROSGP algorithm

    end
   
   
mustar(i,:) = mu;
varstar(i,:) = s2;
m_hist(i,:) = norm(model.postm);
S_hist(i,:) = norm(model.postS);
   
   
   
end



figure;
hold on; grid on;box on;
fillColor = [236 241 249]/255;
fill([Xtest; Xtest(end:-1:1)], [mustar;mustar(end:-1:1)]...
               + 2*[sqrt(varstar); -sqrt(varstar(end:-1:1))], fillColor,'EdgeColor','none');
plot(xm,ym,'k.');
plot(Xtest, mustar, 'b-.', 'LineWidth', 2) % mean predictions in blue
plot(Xtest, ff', 'r', 'LineWidth', 2) % 
plot(Xtest, mustar + 2*sqrt(varstar),'-','LineWidth',1,'color',fillColor) % plus/minus 2 std deviation
plot(Xtest, mustar - 2*sqrt(varstar),'-','LineWidth',1,'color',fillColor) % plus/minus 2 std deviation
plot(model.Xu,-3*ones(size(model.Xu)),'b+','linewidth',2,'markersize',10);
plot(Xuinit,-3*ones(size(Xuinit)),'r+','linewidth',2,'markersize',10);
set(gca, 'fontsize', 18);
box on;grid on;
legend('95% Conf','Training Data','Predictive Mean','True Function');



