clc;clear;

load res_stack.mat

load training_in_p.mat
load training_out_p.mat
load training_target_p.mat    

Max_iter = 25;
test = 1;
flag_init = 1; % 1-> initialize the sGP  the same as traiend lGP except for posterior
%% add new data to the data set
H = 5;
idx = [];
res_2nd_stack = cell(Max_iter,1);

for ii = 1: Max_iter
        ii 

    if isempty(cell2mat(res_stack(ii))) ~= 1

         res = res_stack{ii};
         x   = res.x;
         u   = res.u;
         state_dot_hist = res.state_dot_hist;
         para = res.para;
         init_LTI;

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
              Smodel_y.postS = eye(Nu)*1e2;
        
              Smodel_z.postm = zeros(Smodel_z.m,1);
              Smodel_z.postS = eye(Nu)*1e2;
        
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
            
            logtheta0_z(1:D,1) = 2*ones(D,1); % smoother estimation
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
        save('training_in_1.mat','training_in_p');
        save('training_out_1.mat','training_out_p');
        save('training_target_1.mat','training_target_p');


        try
           res_2 = DGPMPC_After_train_2nd;
               
        catch WARN
            warning('Not converged!');
            idx = [idx,ii];

        continue;
             
         end

    res_2nd_stack{ii,1}  = res_2;



       
   end



 end


%%

x = res_2nd_stack{2}.x;
xd = res.para.path; xd = xd';
H = 5;

green = [167,202,122]/255;
red = [216 56 58]/255;
yellow = [244,205,104]/255;
darkred = [255,122,66]/255;
blue = [47 127 193]/255;
% [95 151 210]/255;


% wind_direction =winds(1,:)/norm(winds(1,:))/2;
figure;
plot3(x(1,1:end-H+1),x(2,1:end-H+1),x(3,1:end-H+1),'color',blue,'linewidth',2)
grid on; hold on; box on 
plot3(xd(1,1:end-H+1),xd(2,1:end-H+1),xd(3,1:end-H+1),'-.','linewidth',2,'color',red);
plot3(x(1,1),x(2,1),x(3,1),'o','color',darkred,'Markersize',20,'MarkerFaceColor',darkred);hold on 
plot3(xd(1,end-H+1),xd(2,end-H+1),xd(3,end-H+1),'p','color',yellow,'Markersize',20,'MarkerFaceColor',yellow);hold on 
plot3(x(1,end-H+1),x(2,end-H+1),x(3,end-H+1),'p','color',green,'Markersize',20,'MarkerFaceColor',green);hold on 
% quiver3(0.2,0.5,3,wind_direction(1),wind_direction(2),wind_direction(3),'linewidth',2,'color','m','MaxHeadSize',0.26);
legend('Actual','Ref','Initial','Final-desired','Final-actual');
xlabel('X/m');
ylabel('Y/m');
zlabel('Z/m');
view(60,10);
set(gca,'fontsize',16,'Fontname','times new roman');
