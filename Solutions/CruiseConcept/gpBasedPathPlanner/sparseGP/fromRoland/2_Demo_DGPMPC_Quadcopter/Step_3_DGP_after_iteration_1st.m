%% train the static GP after each iteration
% When each iteration finishes, we have a set of new data
% add the new data points to the static GP training set and retrain it 
%% GP input : x y 
%  GP output: xdot - f(x,u)-mu_d
%% warning!! 此文档不可连续多次运行--  (运行完此文档，若运行mpc发散，不可再回来重新训练)
% 解决：把1st，2st，3rd 的training in/out  分别保存记录
close all
load training_in_p.mat
load training_out_p.mat
load training_target_p.mat    
% load para.mat
test = 1;
flag_init = 0; % 1-> initialize the sGP  the same as traiend lGP except for posterior
%% add new data to the data set
H = para.trans.H;
training_in_p_new = getTrain_P(x,u,para); % phi,theta,psi,vx,vy,vz,F 
training_out_p_new = getDiff_p(x,u,state_dot_hist,para,LTI);

training_in_p_new = training_in_p_new(150:end,:);
training_out_p_new = training_out_p_new(150:end,:);

[nn,~] = size(training_in_p_new);

% argument the new & old  data together
training_in_1 = [training_in_p;training_in_p_new];
[N,D] = size(training_in_1);
sigma1 = 1e-2*0.268;sigma2 = 1e-2*0.640;sigma3 = 1e-3*0.114;
noise = [sigma1*randn(nn,1) sigma2*randn(nn,1) sigma3*randn(nn,1)];
training_target_p_new =  training_out_p_new+noise;
training_out_1 = [training_out_p; training_out_p_new];
training_target_1 = [training_target_p; training_target_p_new];

figure;
plot(training_target_1);   % 尽量不要有抖振

%% batch train the long term GP
Nu = 20;
options.Likelihood = 'Gaussian';
options.m = Nu; 

% type of inducing variables
options.indType = 'pseudoIns'; % 'pseudoIns' or 'weights'
options.objectFunc = 'var';

new_model_x = varsgpCreate('seard', training_in_1, training_target_1(:,1), options);
new_model_y = varsgpCreate('seard', training_in_1, training_target_1(:,2), options);
new_model_z = varsgpCreate('seard', training_in_1, training_target_1(:,3), options);

index  = randsample(N,new_model_x.m);

Xuinit = training_in_1(index,:);  % this is the initial inducing points, Nu*D


new_model_x.Xu = Xuinit; 
new_model_y.Xu = Xuinit;
new_model_z.Xu = Xuinit;

trops(1) = 500; % number of iterations
trops(2) = 1; % type of optimizations
trops(3) = 1; 

% initialization of the model hyperparameters 
logtheta0_x(1:D,1) = log((max(training_in_1)-min(training_in_1))'/2);
logtheta0_x(D+1,1) = 0.5*log(var(training_out_1(:,1),1)); 
logtheta0_x(D+2,1) = 0.5*log(var(training_out_1(:,1),1)/4);  
new_model_x.GP.logtheta = logtheta0_x(1:end-1);
new_model_x.Likelihood.logtheta = logtheta0_x(end);

logtheta0_y(1:D,1) = log((max(training_in_1)-min(training_in_1))'/2);
logtheta0_y(D+1,1) = 0.5*log(var(training_out_1(:,2),1)); 
logtheta0_y(D+2,1) = 0.5*log(var(training_out_1(:,2),1)/4);  
new_model_y.GP.logtheta = logtheta0_y(1:end-1);
new_model_y.Likelihood.logtheta = logtheta0_y(end);

logtheta0_z(1:D,1) = log((max(training_in_1)-min(training_in_1))'/2);
logtheta0_z(D+1,1) = 0.5*log(var(training_out_1(:,3),1)); 
logtheta0_z(D+2,1) = 0.5*log(var(training_out_1(:,3),1)/4);  
new_model_z.GP.logtheta = logtheta0_z(1:end-1);
new_model_z.Likelihood.logtheta = logtheta0_z(end);

%% GP offline train and predict

[new_model_x, margLogL_x] = varsgpTrain(new_model_x, trops);
[new_model_y, margLogL_y] = varsgpTrain(new_model_y, trops);
[new_model_z, margLogL_z] = varsgpTrain(new_model_z, trops);

for i = 1:N

    [mu_x(i,:), var_x(i,:)] = varsgpPredict_post_x(new_model_x,training_in_1(i,:));
    [mu_y(i,:), var_y(i,:)] = varsgpPredict_post_y(new_model_y,training_in_1(i,:));
    [mu_z(i,:), var_z(i,:)] = varsgpPredict_post_z(new_model_z,training_in_1(i,:));

end

std_x  = sqrt(var_x);
std_y  = sqrt(var_y);
std_z  = sqrt(var_z);
%%
figure;
subplot(3,1,1)
patch([1:N, fliplr(1:N)],[mu_x-2*std_x; flipud(mu_x+2*std_x)],[236 241 249]/255, 'EdgeColor', 'none'); hold on 
plot(mu_x);
plot(training_out_1(:,1));
legend('95% conf.','esti','true');

subplot(3,1,2)
patch([1:N, fliplr(1:N)],[mu_y-2*std_y; flipud(mu_y+2*std_y)],[236 241 249]/255, 'EdgeColor', 'none'); hold on 
plot(mu_y);
plot(training_out_1(:,2));
legend('95% conf.','esti','true');

subplot(3,1,3)
patch([1:N, fliplr(1:N)],[mu_z-2*std_z; flipud(mu_z+2*std_z)],[236 241 249]/255, 'EdgeColor', 'none'); hold on 
plot(mu_z);
plot(training_out_1(:,3));
legend('95% conf.','esti','true');


model_x = new_model_x;
model_y = new_model_y;
model_z = new_model_z;

%% re-initialize the short term GP
if flag_init ==  1

      Smodel_x = model_x;
      Smodel_y = model_y;
      Smodel_z = model_z;

      Smodel_x.postm = zeros(Smodel_x.m,1);
      Smodel_x.postS = eye(Nu)*1e2;

      Smodel_y.postm = zeros(Smodel_y.m,1);
      Smodel_y.postS = eye(Nu)*1e1;

      Smodel_z.postm = zeros(Smodel_z.m,1);
      Smodel_z.postS = eye(Nu)*1e1;

else

    new_Smodel_x = varsgpCreate('seard', training_in_1, training_target_1(:,1), options);
    new_Smodel_y = varsgpCreate('seard', training_in_1, training_target_1(:,2), options);
    new_Smodel_z = varsgpCreate('seard', training_in_1, training_target_1(:,3), options);
    
    index  = randsample(nn,new_Smodel_x.m);  % 换一个更均匀的方式
    
    % Xuinit = training_in_1(index,:);  % this is the initial inducing points, Nu*D
    % Xuinit = training_in_p_new(index,:);  
    
    new_Smodel_x.Xu = Xuinit; 
    new_Smodel_y.Xu = Xuinit;
    new_Smodel_z.Xu = Xuinit;
    
    logtheta0_x(1:D,1) = 2*ones(D,1); % smoother estimation
    logtheta0_x(D+1,1) = -2;
    logtheta0_x(D+2,1) = -5;
    new_Smodel_x.GP.logtheta = logtheta0_x(1:end-1);
    new_Smodel_x.Likelihood.logtheta = logtheta0_x(end);
    
    logtheta0_y(1:D,1) = 2.5*ones(D,1); % smoother estimation
    logtheta0_y(D+1,1) = -2;
    logtheta0_y(D+2,1) = -5;
    new_Smodel_y.GP.logtheta = logtheta0_y(1:end-1);
    new_Smodel_y.Likelihood.logtheta = logtheta0_y(end);
    
    logtheta0_z(1:D,1) = 1*ones(D,1); % smoother estimation
    logtheta0_z(D+1,1) = -2;
    logtheta0_z(D+2,1) = -5;
    new_Smodel_z.GP.logtheta = logtheta0_z(1:end-1);
    new_Smodel_z.Likelihood.logtheta = logtheta0_z(end);
    
    new_Smodel_x.postm = zeros(new_Smodel_x.m,1);
    new_Smodel_x.postS = eye(Nu)*1e2;
    [~,~,Knm,invKm,L,Kmm] = computeKminv(new_Smodel_x, training_in_1);
    new_Smodel_x.Kmm = Kmm;
    new_Smodel_x.Kmn = Knm';
    new_Smodel_x.Knm = Knm;
    new_Smodel_x.L   = L;
    new_Smodel_x.invKm = invKm;
    
    new_Smodel_y.postm = zeros(new_Smodel_y.m,1);
    new_Smodel_y.postS = eye(Nu)*1e2;
    [~,~,Knm,invKm,L,Kmm] = computeKminv(new_Smodel_y, training_in_1);
    new_Smodel_y.Kmm = Kmm;
    new_Smodel_y.Kmn = Knm';
    new_Smodel_y.Knm = Knm;
    new_Smodel_y.L   = L;
    new_Smodel_y.invKm = invKm;
    
    new_Smodel_z.postm = zeros(new_Smodel_z.m,1);
    new_Smodel_z.postS = eye(Nu)*1e2;
    [~,~,Knm,invKm,L,Kmm] = computeKminv(new_Smodel_z, training_in_1);
    new_Smodel_z.Kmm = Kmm;
    new_Smodel_z.Kmn = Knm';
    new_Smodel_z.Knm = Knm;
    new_Smodel_z.L   = L;
    new_Smodel_z.invKm = invKm;
    
    
    Smodel_x = new_Smodel_x;
    Smodel_y = new_Smodel_y;
    Smodel_z = new_Smodel_z;  
end


%% resave all the data
training_in_p = training_in_1;
training_out_p= training_out_1;
training_target_p = training_target_1;

save('gpmodel_init.mat','model_x','model_y','model_z','Smodel_x','Smodel_y','Smodel_z');


%% test short term recursion
if test == 1

for i = 1:N

    [mu_x, var_x] = varsgpPredict_post_1(Smodel_x,training_in_1(i,:));
    [mu_y, var_y] = varsgpPredict_post_1(Smodel_y,training_in_1(i,:));
    [mu_z, var_z] = varsgpPredict_post_1(Smodel_z,training_in_1(i,:));

   err_x = training_out_1(i,1)-mu_x;
   err_y = training_out_1(i,2)-mu_y;
   err_z = training_out_1(i,3)-mu_z;

   [phix,~,~,~,~] = computeKminv(Smodel_x, training_in_1(i,:));  % new kernel
   [phiy,~,~,~,~] = computeKminv(Smodel_y, training_in_1(i,:));
   [phiz,~,~,~,~] = computeKminv(Smodel_z, training_in_1(i,:));
        lambda = 0.995;
       Smodel_x = ROSV_GP_lambda(Smodel_x,err_x,phix,lambda);
       Smodel_y = ROSV_GP_lambda(Smodel_y,err_y,phiy,lambda);
       Smodel_z = ROSV_GP_lambda(Smodel_z,err_z,phiz,lambda);



        mustar(i,:) = [mu_x mu_y mu_z];
        varstar(i,:) = [var_x var_y var_z];
%          m_hist(i,:) = norm(Smodel_x.postm);
%         S_hist(i,:) = trace(Smodel_x.postS);

end
 
% varstar(1:10,1:10) = zeros(10,10);


std_x = sqrt(varstar(:,1)); std_y = sqrt(varstar(:,2)); std_z = sqrt(varstar(:,3));


figure;
subplot(3,1,1); hold on 
plot(mustar(:,1));
plot(training_out_1(:,1));
legend('95% conf.','esti','true');

subplot(3,1,2)
hold on 
plot(mustar(:,2));
plot(training_out_1(:,2));
legend('95% conf.','esti','true');

subplot(3,1,3)
hold on 
plot(mustar(:,3));
plot(training_out_1(:,3));
legend('95% conf.','esti','true');


Err_x = sum((mustar(:,1)-training_out_1(:,1)).^2);
MSE_x = Err_x/N;

Err_y = sum((mustar(:,2)-training_out_1(:,2)).^2);
MSE_y = Err_y/N;

Err_z = sum((mustar(:,3)-training_out_1(:,3)).^2);
MSE_z = Err_z/N;

MSE = [MSE_x MSE_y MSE_z]

  
end


save('training_in_1.mat','training_in_p');
save('training_out_1.mat','training_out_p');
save('training_target_1.mat','training_target_p');

%% regenerate DGP mex
codegen DGP_varsgpPredict_x_2 -args {model_x,Smodel_x,zeros(1,7),zeros(7)}

codegen DGP_varsgpPredict_y_2 -args {model_y,Smodel_y,zeros(1,7),zeros(7)}

codegen DGP_varsgpPredict_z_2 -args {model_z,Smodel_z,zeros(1,7),zeros(7)}

