%% An example for the dictionary based online GP

warning off
clc;clear; close all

%  global noise sigma


addpath('misc');
% addpath('SPGP_dist');
addpath('docs');

flag_update = 1;
if_2nd_round = 0;

dT = 0.05;
% dT = 0.02;

minX = -5;
maxX = 0;           % training set
minX_Test = -5;
maxX_Test = 10;      % testing set

fx = @(x) 2*sin(x)+0.05*x.^2+2;
fg = @(x) fx(x)-2*cos(x)-0.2*x;
%  fg = @(x) fx(x)-0.5*cos(x)-0.1*x.^2-0.2*x;

xn = minX_Test:dT:maxX_Test;  % column
N = length(xn);
xn1 = xn(1:ceil(N/2)); xn2 = xn(ceil(N/2)+1:end);
f1 = fx(xn1);
f0 = fx(xn);
f2 = fg(xn2);

ff = [f1 f2];



% choose the measurement data among xn as the training set

Nm = 500;    % training set size  （big）
% xm = sort(minX + rand(1,Nm)*(maxX-minX));  % choose randmoly
% % ym = fx(xm)+ sigma*randn(1,Nm);  % noisy data 
% % Xtrain = xm';
% Ytrain = ym';
xn = minX_Test:dT:maxX_Test;  % column'
Xtest = xn';



%% dic GP  parameters
bandwidth = 0.1;
noise = sqrt(0.5);
% noise = 0.02;
tol =  0.001;  % larger: more tolerence
max_points = 200; % max size of GP dictionary



gpr = onlineGP(bandwidth,noise,max_points,tol);



s22 = tic;
%% sim start
%% 1st round
for i  = 1:N
     

    if i == 1
     
      % initialize gp model
  
      gpr.process(Xtest(i),ff(i));
      [mu,s2] = gpr.predict(Xtest(i));

%     gpr.process(Xu',f_xu);

    else
         [mu,s2] = gpr.predict(Xtest(i));
    end
   
    
     if  flag_update == 1
         t3 = tic;
      gpr.update(Xtest(i), ff(i));
      s3(i)=toc(t3);
   
    end
 
   
mustar(i,:) = mu;
varstar(i,:) = s2;

  alpha = gpr.get('alpha');
   
end

BV = gpr.get('BV');


figure;
plot(s3);
