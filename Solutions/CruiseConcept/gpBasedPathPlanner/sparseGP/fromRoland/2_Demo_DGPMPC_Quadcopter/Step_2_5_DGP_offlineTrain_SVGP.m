%% To learn the difference between nominal model and full model using GP.
% Using modified varGP toolbox

close all; clear; clc;
% set_parameters;
% init_LTI;
clear gpmodel_trans

load training_data.mat
test = 1;  % if test the initial sGP revursive performance
flag_init =0 ; % 1-> initialize the sGP  the same as traiend lGP except for posterior
Nu = 20;  % the number of pseudo inputs
% Nu = 50;
%   Nu = 70;
% Nu = 100;

%% Translational dynamics
training_in_p = getTrain_P(x,u,para); % phi,theta,psi,vx,vy,vz,F 
[N,D] = size(training_in_p);

% sigma1 = 1e-1*0.268;
% sigma2 = 1e-4*0.640;
% sigma3 = 1e-1*0.114;

sigma1 = 1e-2*0.268;
sigma2 = 1e-2*0.640;
sigma3 = 1e-3*0.114;

noise = [sigma1*randn(N,1) sigma2*randn(N,1) sigma3*randn(N,1)];

training_out_p = getDiff_p(x,u,state_dot_hist,para,LTI);
training_target_p =  training_out_p+noise;

%% initialize GP 
options.Likelihood = 'Gaussian';
options.m = Nu; 

% type of inducing variables
options.indType = 'pseudoIns'; % 'pseudoIns' or 'weights'
options.objectFunc = 'var';

model_x = varsgpCreate('seard', training_in_p, training_target_p(:,1), options);
model_y = varsgpCreate('seard', training_in_p, training_target_p(:,2), options);
model_z = varsgpCreate('seard', training_in_p, training_target_p(:,3), options);

index  = randsample(N,model_x.m);

Xuinit = training_in_p(index,:);  % this is the initial inducing points, Nu*D


model_x.Xu = Xuinit; 
model_y.Xu = Xuinit;
model_z.Xu = Xuinit;

trops(1) = 200; % number of iterations
trops(2) = 1; % type of optimizations
trops(3) = 1; 

% initialization of the model hyperparameters 
logtheta0_x(1:D,1) = log((max(training_in_p)-min(training_in_p))'/2);
logtheta0_x(D+1,1) = 0.5*log(var(training_out_p(:,1),1)); 
logtheta0_x(D+2,1) = 0.5*log(var(training_out_p(:,1),1)/4);  
model_x.GP.logtheta = logtheta0_x(1:end-1);
model_x.Likelihood.logtheta = logtheta0_x(end);

logtheta0_y(1:D,1) = log((max(training_in_p)-min(training_in_p))'/2);
logtheta0_y(D+1,1) = 0.5*log(var(training_out_p(:,2),1)); 
logtheta0_y(D+2,1) = 0.5*log(var(training_out_p(:,2),1)/4);  
model_y.GP.logtheta = logtheta0_y(1:end-1);
model_y.Likelihood.logtheta = logtheta0_y(end);

logtheta0_z(1:D,1) = log((max(training_in_p)-min(training_in_p))'/2);
logtheta0_z(D+1,1) = 0.5*log(var(training_out_p(:,3),1)); 
logtheta0_z(D+2,1) = 0.5*log(var(training_out_p(:,3),1)/4);  
model_z.GP.logtheta = logtheta0_z(1:end-1);
model_z.Likelihood.logtheta = logtheta0_z(end);

% Smodel_x = model_x;
% Smodel_y = model_y;
% Smodel_z = model_z;
%% GP offline train and predict

[model_x, margLogL_x] = varsgpTrain(model_x, trops);
[model_y, margLogL_y] = varsgpTrain(model_y, trops);
[model_z, margLogL_z] = varsgpTrain(model_z, trops);

for i = 1:N

    [mu_x(i,:), var_x(i,:)] = varsgpPredict_post_x(model_x,training_in_p(i,:));
    [mu_y(i,:), var_y(i,:)] = varsgpPredict_post_y(model_y,training_in_p(i,:));
    [mu_z(i,:), var_z(i,:)] = varsgpPredict_post_z(model_z,training_in_p(i,:));

end

std_x  = sqrt(var_x);
std_y  = sqrt(var_y);
std_z  = sqrt(var_z);
%%
figure;
subplot(3,1,1)
patch([1:N, fliplr(1:N)],[mu_x-2*std_x; flipud(mu_x+2*std_x)],[236 241 249]/255, 'EdgeColor', 'none'); hold on 
plot(mu_x);
plot(training_out_p(:,1));
legend('95% conf.','esti','true');

subplot(3,1,2)
patch([1:N, fliplr(1:N)],[mu_y-2*std_y; flipud(mu_y+2*std_y)],[236 241 249]/255, 'EdgeColor', 'none'); hold on 
plot(mu_y);
plot(training_out_p(:,2));
legend('95% conf.','esti','true');

subplot(3,1,3)
patch([1:N, fliplr(1:N)],[mu_z-2*std_z; flipud(mu_z+2*std_z)],[236 241 249]/255, 'EdgeColor', 'none'); hold on 
plot(mu_z);
plot(training_out_p(:,3));
legend('95% conf.','esti','true');


% model_x.postS = eye(Nu)*1e0;   model_y.postS = eye(Nu)*1e0;    model_z.postS = eye(Nu)*1e0;
% model_x.postS0 = model_x.postS; model_y.postS0 = model_y.postS; model_z.postS0 = model_z.postS;

% 
% Smodel_x = model_x;
% Smodel_y = model_y;
% Smodel_z = model_z;% initializes as untrained GP

%% initialize the short term gp

if flag_init ==  1

      Smodel_x = model_x;
      Smodel_y = model_y;
      Smodel_z = model_z;

      Smodel_x.postm = zeros(Smodel_x.m,1);
      Smodel_x.postS = eye(Nu)*1e2;

      Smodel_y.postm = zeros(Smodel_y.m,1);
      Smodel_y.postS = eye(Nu)*1e3;

      Smodel_z.postm = zeros(Smodel_z.m,1);
      Smodel_z.postS = eye(Nu)*1e2;

      fields = {'alpha','sigma2'};

      Smodel_x = rmfield(Smodel_x,fields);
      Smodel_y = rmfield(Smodel_y,fields);
      Smodel_z = rmfield(Smodel_z,fields);

else

        Smodel_x = varsgpCreate('seard', training_in_p, training_target_p(:,1), options);
        Smodel_y = varsgpCreate('seard', training_in_p, training_target_p(:,2), options);
        Smodel_z = varsgpCreate('seard', training_in_p, training_target_p(:,3), options);
        
        index  = randsample(N,Smodel_x.m);  % 换一个更均匀的方式
        
        % Xuinit = mean(training_in_p)
        
        Xuinit = training_in_p(index,:);  % this is the initial inducing points, Nu*D
        
        Smodel_x.Xu = Xuinit; 
        Smodel_y.Xu = Xuinit;
        Smodel_z.Xu = Xuinit;
        
        logtheta0_x(1:D,1) = 2*ones(D,1); % smoother estimation
        logtheta0_x(D+1,1) = -2;
        logtheta0_x(D+2,1) = -5;
        Smodel_x.GP.logtheta = logtheta0_x(1:end-1);
        Smodel_x.Likelihood.logtheta = logtheta0_x(end);
        
        logtheta0_y(1:D,1) = 2.5*ones(D,1); % smoother estimation
        logtheta0_y(D+1,1) = -2;
        logtheta0_y(D+2,1) = -5;
        Smodel_y.GP.logtheta = logtheta0_y(1:end-1);
        Smodel_y.Likelihood.logtheta = logtheta0_y(end);
        
        logtheta0_z(1:D,1) =2*ones(D,1); % smoother estimation
        logtheta0_z(D+1,1) = -2;
        logtheta0_z(D+2,1) = -5;
        Smodel_z.GP.logtheta = logtheta0_z(1:end-1);
        Smodel_z.Likelihood.logtheta = logtheta0_z(end);
        
        Smodel_x.postm = zeros(Smodel_x.m,1);
        Smodel_x.postS = eye(Nu)*1e2;
        [~,~,Knm,invKm,L,Kmm] = computeKminv(Smodel_x, training_in_p);
        Smodel_x.Kmm = Kmm;
        Smodel_x.Kmn = Knm';
        Smodel_x.Knm = Knm;
        Smodel_x.L   = L;
        Smodel_x.invKm = invKm;
        
        Smodel_y.postm = zeros(Smodel_y.m,1);
        Smodel_y.postS = eye(Nu)*1e2;
        [~,~,Knm,invKm,L,Kmm] = computeKminv(Smodel_y, training_in_p);
        Smodel_y.Kmm = Kmm;
        Smodel_y.Kmn = Knm';
        Smodel_y.Knm = Knm;
        Smodel_y.L   = L;
        Smodel_y.invKm = invKm;
        
        Smodel_z.postm = zeros(Smodel_z.m,1);
        Smodel_z.postS = eye(Nu)*1e2;
        [~,~,Knm,invKm,L,Kmm] = computeKminv(Smodel_z, training_in_p);
        Smodel_z.Kmm = Kmm;
        Smodel_z.Kmn = Knm';
        Smodel_z.Knm = Knm;
        Smodel_z.L   = L;
        Smodel_z.invKm = invKm;
end


save('gpmodel_init.mat','model_x','model_y','model_z','Smodel_x','Smodel_y','Smodel_z');
save('training_in_p.mat','training_in_p');
save('training_out_p.mat','training_out_p');
save('training_target_p.mat','training_target_p');

%% test short term recursion
if test == 1

for i = 1:N

    [mu_x, var_x] = varsgpPredict_post_1(Smodel_x,training_in_p(i,:));
    [mu_y, var_y] = varsgpPredict_post_1(Smodel_y,training_in_p(i,:));
    [mu_z, var_z] = varsgpPredict_post_1(Smodel_z,training_in_p(i,:));

   err_x = training_out_p(i,1)-mu_x;
   err_y = training_out_p(i,2)-mu_y;
   err_z = training_out_p(i,3)-mu_z;

   lambda = 0.99;
   [phix,~,~,~,~] = computeKminv(Smodel_x, training_in_p(i,:));  % new kernel
   [phiy,~,~,~,~] = computeKminv(Smodel_y, training_in_p(i,:));
   [phiz,~,~,~,~] = computeKminv(Smodel_z, training_in_p(i,:));
         
       Smodel_x = ROSV_GP_lambda(Smodel_x,err_x,phix,lambda);
       Smodel_y = ROSV_GP_lambda(Smodel_y,err_y,phiy,lambda);
       Smodel_z = ROSV_GP_lambda(Smodel_z,err_z,phiz,lambda);

% 
%         S = Smodel_x.postS;
%         m = Smodel_x.postm;
% 
%         Gk = lambda + phi*S*phi';
% 
%         iGk = 1/Gk;
% 
%         Lk   = iGk * S * phi';   % 20*1
% 
%         Hk = S * phi' * iGk;
% 
%         S_new = lambda^(-1)*(S -  Hk*Gk*Hk');
% 
% 
%         Smodel_x.postS = S_new;
%         Smodel_x.postm = m + Lk*err;

        mustar(i,:) = [mu_x mu_y mu_z];
        varstar(i,:) = [var_x var_y var_z];
%          m_hist(i,:) = norm(Smodel_x.postm);
%         S_hist(i,:) = trace(Smodel_x.postS);

end
 
varstar(1:10,1:10) = zeros(10,10);


std_x = sqrt(varstar(:,1)); std_y = sqrt(varstar(:,2)); std_z = sqrt(varstar(:,3));


figure;
subplot(3,1,1); hold on 
% patch([1:N, fliplr(1:N)],[mustar(:,1)-2*std_x; flipud(mustar(:,1)+2*std_x)],[236 241 249]/255, 'EdgeColor', 'none');
plot(mustar(:,1));
plot(training_out_p(:,1));
legend('95% conf.','esti','true');

subplot(3,1,2)
% patch([1:N, fliplr(1:N)],[mustar(:,2)-2*std_y; flipud(mustar(:,2)+2*std_y)],[236 241 249]/255, 'EdgeColor', 'none'); 
hold on 
plot(mustar(:,2));
plot(training_out_p(:,2));
legend('95% conf.','esti','true');

subplot(3,1,3)
% patch([1:N, fliplr(1:N)],[mustar(:,3)-2*std_z; flipud(mustar(:,3)+2*std_z)],[236 241 249]/255, 'EdgeColor', 'none');
hold on 
plot(mustar(:,3));
plot(training_out_p(:,3));
legend('95% conf.','esti','true');


Err_x = sum((mustar(:,1)-training_out_p(:,1)).^2);
MSE_x = Err_x/N;

Err_y = sum((mustar(:,2)-training_out_p(:,2)).^2);
MSE_y = Err_y/N;

Err_z = sum((mustar(:,3)-training_out_p(:,3)).^2);
MSE_z = Err_z/N;

MSE = [MSE_x MSE_y MSE_z]

  
end





