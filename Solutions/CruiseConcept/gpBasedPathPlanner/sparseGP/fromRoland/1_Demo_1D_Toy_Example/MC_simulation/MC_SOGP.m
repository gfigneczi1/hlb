%% An example for the dictionary based online GP
function s3 = MC_SOGP(tol,bandwidth)
warning off
clc; close all

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
% bandwidth = 2;
noise = sqrt(0.5);
% noise = 0.02;
% tol =  0.001;  % larger: more tolerence
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

%% second round
if if_2nd_round == 1
for i  = 1:N
     
% in this round, no initialization. use the BV from last round
    
         [mu,s2] = gpr.predict(Xtest(i));
   
   
    
     if  flag_update == 1
         
      gpr.update(Xtest(i), ff(i));
   
    end
 
   
mustar_2nd(i,:) = mu;
varstar_2nd(i,:) = s2;

  alpha = gpr.get('alpha');
   
end

end

BV = gpr.get('BV');


t2= toc(s22)

% figure;
% plot(xn,mustar,'b-','linewidth',2);hold on
% plot(xn,ff','r:','linewidth',2);
% plot(BV,zeros(1,length(BV)),'r+','linewidth',2);hold on
% legend('1st round');
% 
% figure;
% plot(xn,mustar_2nd,'b-','linewidth',2);hold on
% plot(xn,ff','r:','linewidth',2);
% plot(BV,zeros(1,length(BV)),'r+','linewidth',2);hold on
% legend('2nd round');
% 
% % figure;
% % plot(varstar);
% % legend('Predicted variance');
% 

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

% plot(model.Xu,-3*ones(size(model.Xu)),'b+','linewidth',2,'markersize',10);
% plot(Xuinit,-3*ones(size(Xuinit)),'r+','linewidth',2,'markersize',10);
% axis([min(Xtest) max(Xtest) -3 3]);
set(gca, 'fontsize', 18);
%set(gca, 'YTick', []);
%set(gca, 'XTick', []);
box on;grid on;
legend('95% Conf','Training data','Mean','True');

if if_2nd_round == 1
figure;
hold on;
fillColor = [236 241 249]/255;
fill([Xtest; Xtest(end:-1:1)], [mustar_2nd;mustar_2nd(end:-1:1)]...
               + 2*[sqrt(varstar_2nd); -sqrt(varstar_2nd(end:-1:1))], fillColor,'EdgeColor','none');
hold on; grid on;box on;
% plot(Xtest,muFull,'linewidth',1.5);
%  plot(Xtrain,Ytrain,'k.','linewidth',1);
% plot(Xtest,fx(Xtest)-my,'linewidth',1.5);
set(gca,'Fontsize',15);

      
% plot(Xtrain, Ytrain,'.y','markersize',5); % data points in magenta

plot(Xtest, mustar_2nd, 'b-.', 'LineWidth', 2) % mean predictions in blue
plot(Xtest, ff', 'r', 'LineWidth', 2) % 
plot(Xtest, mustar_2nd + 2*sqrt(varstar_2nd),'-','LineWidth',1,'color',fillColor) % plus/minus 2 std deviation
plot(Xtest, mustar_2nd - 2*sqrt(varstar_2nd),'-','LineWidth',1,'color',fillColor) % plus/minus 2 std deviation

% plot(model.Xu,-3*ones(size(model.Xu)),'b+','linewidth',2,'markersize',10);
% plot(Xuinit,-3*ones(size(Xuinit)),'r+','linewidth',2,'markersize',10);
% axis([min(Xtest) max(Xtest) -3 3]);
set(gca, 'fontsize', 18);
%set(gca, 'YTick', []);
%set(gca, 'XTick', []);
box on;grid on;
legend('95% Conf','Training data','Mean','True');


save('dic_GP.mat','mustar','mustar_2nd','varstar','varstar_2nd');

end