%% lGPMPC- with online GP compensation in prediction model
% continous time-varying wind added in dynamics
% xdot = Ax+Bu , nominal model
% collect data using the random trajectory
% only with the long-term GP
%%
clc;clear all;close all
global winds flag_wind  i

warning off
addpath('util');
addpath('odefiles');
addpath('varGP');
addpath('GP_related');
addpath('MPCsolver');

% load gpmodel_trans.mat
load gpmodel_init.mat

set_parameters;
init_LTI;

if_gp_on = 1;
flag_update = 0;
flag_wind = 1; % if there is wind disturbances or not
winds = generate_wind(para,3)*flag_wind;

tt = 0:0.05:20;
figure;
plot(tt',winds(:,1),'r:','linewidth',2); grid on; hold on 
plot(tt',winds(:,2),'b-','linewidth',2);
plot(tt',winds(:,3),'k--','linewidth',2);
legend('x','y','z');
xlabel('Times/s');
ylabel('Wind Speed m/s');
set(gca,'fontsize',16,'Fontname','times new roman');



%% Dyn settings

para.trans.A  = LTI_trans.A;
para.trans.B  = LTI_trans.B;  % discrete 

para.rot.A  =  LTI_rot.A;
para.rot.B  =  LTI_rot.B;  % discrete 

para.model_x = model_x;
para.model_y = model_y;
para.model_z = model_z;

% para.Smodel_x = Smodel_x;
% para.Smodel_y = Smodel_y;
% para.Smodel_z = Smodel_z;


[~,K_rot,~] = idare(para.rot.A,para.rot.B,para.rot.Q,para.rot.R);
kp = diag([2,2,2]);
kd = diag([0.5,0.5,0.5]);

% attitude control
kp_roll_pitch_yaw = [6.0, 6.0,4.0];
kd_roll_pitch_yaw = [8.0 8.0 8.0];
%%
ref = generate_Ref(para,2);  % 1--fixed point 2--helix 3--''8'' 4--random
xd = ref.x(1:6,:);
para.path  = xd';

xtrans0 = xd(:,1);
% xtrans0 = zeros(6,1);
utrans0 = zeros(3,1);
xrot0 = zeros(6,1);
urot0 = zeros(3,1);

N = para.N;
time = para.time;
dT = para.dT;

xtrans = zeros(6,N);
utrans = zeros(3,N);
xtrans(:,1) = xtrans0;

xrot = zeros(6,N);
urot = zeros(3,N);
xrot(:,1) = xrot0;
x = [xtrans;xrot];

mu_d_hist  = [];
delta_hist = [];
mu_a = zeros(3,1);
err_hist = [];

tic;
     for i = 1: para.N-para.trans.H

        %-----------------------------------------------------------%
        %   Compute the optimal control with in the horizon         %
        %-----------------------------------------------------------%
          
           
          [u_opt_trans,cost_trans(i)]  = mpcCtrl_trans_lGP(xtrans(:,i),xrot(:,i),time(i),para,i); % considering the variance in prediction model

         
%             [u_opt_trans,cost_trans(i)]  = mpcCtrl_trans_GP_var(xtrans(:,i),time(i),para,i);
          utrans(:,i) =  u_opt_trans;
          
           u1(i)  = utrans(3,i) + para.m*para.g;   % thrust force
           ax_cmd = -utrans(2,i)*para.g;
           ay_cmd = utrans(1,i)*para.g;
           az_cmd = -u1(i)/para.m;   
           
          % retrieve the attitude information from the updated state
             phi   = xrot(1,i);
             theta = xrot(2,i);
             psi   = xrot(3,i);
             p     = xrot(4,i);
             q     = xrot(5,i);
             r     = xrot(6,i);
             
             
           % roll, pitch control
            [p_cmd, q_cmd] = RollPitch_control(u1(i), ax_cmd, ay_cmd, phi, theta, psi, kp_roll_pitch_yaw,para);

           % yaw control
              r_cmd     = Yaw_control(ref.psi(i), psi, kp_roll_pitch_yaw);
                 
           % pqr (body-rate) controller
        [pd_cmd, qd_cmd, rd_cmd] = PQR_control(p_cmd, q_cmd, r_cmd, p, q, r, kd_roll_pitch_yaw);        
                 
        pd_cmd = pd_cmd -  if_gp_on*mu_a(1);
        qd_cmd = qd_cmd -  if_gp_on*mu_a(2);
        rd_cmd = rd_cmd -  if_gp_on*mu_a(3);
      %% Set propeller speeds
          % Set angular velocities on each propeller
           J = [para.Jx para.Jy para.Jz];
         [w1, w2, w3, w4] = setAngularVelocities(u1(i), pd_cmd, qd_cmd, rd_cmd, x(:,i),J,para);
      
      %% Compute the updated state
             f1 = para.Ct*w1^2;   f2 = para.Ct*w2^2;   f3 = para.Ct*w3^2;   f4 = para.Ct*w4^2;
             Ftotal = f1 + f2 + f3 + f4;
        
              F_hist(i) = Ftotal;
     
      % Torques/Moments
             t1 =  para.Cq*w1^2;  t2 = -para.Cq*w2^2;  t3 =  para.Cq*w3^2;  t4 = -para.Cq*w4^2;

             tx = (f1 + f4 - f2 - f3)*para.l;
             ty = (f1 + f2 - f3 - f4)*para.l;
             tz = t1 + t2 + t3 + t4;
             tau = [tx ty tz];
             urot(:,i)  = tau';

        %-----------------------------------------------------------%
        %                Simulate the real model                    %
        %-----------------------------------------------------------%
         
       
        u(:,i) = [utrans(:,i);urot(:,i)];
        
        % data collection : the measurement of state dot
        state_dot = dynsys_Plin_Anonlin_wind(i,x(:,i),u(:,i),LTI,para,i);
        state_dot_hist(:,i) = state_dot;
      
        xnew = ode4(@dynsys_Plin_Anonlin_wind,[i-1,i]*dT,x(:,i),u(:,i),LTI,para,i);
        x(:,i+1) = xnew(end,:)';
        xtrans(:,i+1) =  x(1:6,i+1);
        xrot(:,i+1)   =  x(7:12,i+1);
    
       

        %------------------------------------------------------------%
        %                 Compute the disturbances                   %
        %------------------------------------------------------------%
%          
          mu_xu = [xrot(1:3,i+1);x(4:6,i+1);u(3,i)]';
%         
%         % long term 

        [mu_lp1, var_lp1] = varsgpPredict_post_x(model_x,mu_xu);
        [mu_lp2, var_lp2] = varsgpPredict_post_y(model_y,mu_xu);
        [mu_lp3, var_lp3] = varsgpPredict_post_z(model_z,mu_xu);
%        
       % short term
  
%        [mu_sp1, var_sp1] = varsgpPredict_post_x(Smodel_x,mu_xu);
%        [mu_sp2, var_sp2] = varsgpPredict_post_y(Smodel_y,mu_xu);
%        [mu_sp3, var_sp3] = varsgpPredict_post_z(Smodel_z,mu_xu);



          var_sd = [var_lp1,var_lp2,var_lp3];

          mu_a = zeros(3,1);  %
          mu_lp = [mu_lp1 mu_lp2 mu_lp3]; %1x3
       
          
%           mu_p = mu_lp + mu_sp;
          
          mu_p_hist(i,:) = mu_lp;
          
          mu_lp_hist(i,:) = mu_lp;
%           mu_sp_hist(i,:) = mu_sp;
           var_sd_hist(i,:) = var_sd;
          
%           mu_a_hist(i,:) = mu_a';

         target_p = getDiffOnline_p(x(:,i+1),u(:,i),state_dot_hist(:,i),para,LTI); 
         
%         target_sp  = getDiffOnline_p_s(x(:,i+1),u(:,i),mu_lp,state_dot_hist(:,i),para,LTI);
%          target_sp  = target_p - mu_lp;
         
         target_p_hist(i,:)  = target_p; % 1x3
%          target_sp_hist(i,:) =  target_sp;
       
%          
%          if flag_update == 1
%              
%          [Smodel_x,Smodel_y,Smodel_z,err] = GP_update_SVGP(mu_xu,target_sp,mu_sp,Smodel_x,Smodel_y,Smodel_z); % err 1x3 vector
%          
%          para.Smodel_x = Smodel_x;
%          para.Smodel_y = Smodel_y;
%          para.Smodel_z = Smodel_z;
%          
%          err_hist = [err_hist;err];
%             
%          end

     end
 toc       
%%
H = para.trans.H;

figure;
subplot(3,1,1)
plot(para.time(1,1:end-H+1),x(1,1:end-H+1),'b','linewidth',1.5);
hold on
plot(para.time(1,1:end-H+1),xd(1,1:end-H+1),'r--','linewidth',1.5);
legend('x','xref');

subplot(3,1,2)
plot(para.time(1,1:end-H+1),x(2,1:end-H+1),'b','linewidth',1.5);
hold on
plot(para.time(1,1:end-H+1),xd(2,1:end-H+1),'r--','linewidth',1.5);
legend('y','yref');

subplot(3,1,3)
plot(para.time(1,1:end-H+1),x(3,1:end-H+1),'b','linewidth',1.5);
hold on
plot(para.time(1,1:end-H+1),xd(3,1:end-H+1),'r--','linewidth',1.5);
legend('z','zref');
xlabel('Time/s');


figure;
subplot(3,1,1)
plot(para.time(1,1:end-H+1),x(4,1:end-H+1),'b','linewidth',1.5);
hold on
plot(para.time(1,1:end-H+1),xd(4,1:end-H+1),'r--','linewidth',1.5);
legend('vx','vxref');

subplot(3,1,2)
plot(para.time(1,1:end-H+1),x(5,1:end-H+1),'b','linewidth',1.5);
hold on
plot(para.time(1,1:end-H+1),xd(5,1:end-H+1),'r--','linewidth',1.5);
legend('vy','vyref');

subplot(3,1,3)
plot(para.time(1,1:end-H+1),x(6,1:end-H+1),'b','linewidth',1.5);
hold on
plot(para.time(1,1:end-H+1),xd(6,1:end-H+1),'r--','linewidth',1.5);
legend('vz','vzref');
xlabel('Time/s');




wind_direction =winds(1,:)/norm(winds(1,:))/2;
figure;
plot3(x(1,1:end-H+1),x(2,1:end-H+1),x(3,1:end-H+1),'b','linewidth',1.5)
grid on; hold on; box on 
plot3(xd(1,1:end-H+1),xd(2,1:end-H+1),xd(3,1:end-H+1),'r--','linewidth',1.5);
plot3(x(1,1),x(2,1),x(3,1),'ro','Markersize',20,'MarkerFaceColor',[.49 1 .63]);hold on 
plot3(xd(1,end-H+1),xd(2,end-H+1),xd(3,end-H+1),'p','Markersize',20,'MarkerFaceColor',[.49 1 .63]);hold on 
plot3(x(1,end-H+1),x(2,end-H+1),x(3,end-H+1),'p','Markersize',20,'MarkerFaceColor',[1 0 0]);hold on 
quiver3(0.2,0.5,3,wind_direction(1),wind_direction(2),wind_direction(3),'linewidth',2,'color','m','MaxHeadSize',0.26);
legend('Actual','Ref','Initial','Final-desired','Final-actual','wind');
xlabel('X,m');
ylabel('Y,m');
zlabel('Z,m');
%     title('Flight Trajectory During Data Collection Phase')
view(60,10);
set(gca,'Fontsize',15);


figure;
plot(x(1,1:end-1),x(2,1:end-1),'b-','linewidth',1.5);hold on 
plot(xd(1,1:end-1),xd(2,1:end-1),'r--','linewidth',1.5);
plot(x(1,1),x(2,1),'go','markersize',20,'linewidth',2);
xlabel('x'); ylabel('y');




figure;
plot(x(1,1:30),x(2,1:30),'b-','linewidth',1.5);hold on 
plot(xd(1,1:end-1),xd(2,1:end-1),'r--','linewidth',1.5);
plot(x(1,1),x(2,1),'go','markersize',20,'linewidth',2);
xlabel('x'); ylabel('y');


figure;
plot(time(1:end-H),u(1:3,:),'b-','linewidth',1.5); hold on
legend('Control input-translational');
grid on 
xlabel('Time/s');
set(gca,'Fontsize',16);

figure;
plot(F_hist,'linewidth',1.5);hold on
plot(u1,'linewidth',1.5);
legend('Thrust force')

figure;
plot(time(1:end-H),u(4:6,:),'b-','linewidth',1.5); hold on;
legend('Control input-rotational')
grid on 
xlabel('Time/s');
set(gca,'Fontsize',16);


figure;
plot(time(1:end-H),cost_trans,'b-','linewidth',1.5); hold on ;box on; grid on 
legend('Cost Value');
xlabel('Time/s');
set(gca,'Fontsize',16);


Err_x = sum((x(1,1:end-H)-xd(1,1:end-H)).^2);
MSE_x = Err_x/(N-H);

Err_y = sum((x(2,1:end-H)-xd(2,1:end-H)).^2);
MSE_y = Err_y/(N-H);

Err_z = sum((x(3,1:end-H)-xd(3,1:end-H)).^2);
MSE_z = Err_z/(N-H);

MSE = [MSE_x MSE_y MSE_z]

  
figure('Name','GP estimation result-long');plot(mu_lp_hist);
       
figure('Name','GP estimation result-total ')
subplot(3,1,1)
hold on ;grid on;box on
plot(time(1:end-H),target_p_hist(:,1),'r--','linewidth',1.5);
hold on 
plot(time(1:end-H),mu_p_hist(:,1),'b-','linewidth',1.5);
legend('True','Estimate')
ylabel('$\Delta_x$','interpreter','latex');
set(gca,'fontsize',16,'Fontname','times new roman');


subplot(3,1,2)
hold on ;grid on;box on
plot(time(1:end-H),target_p_hist(:,2),'r--','linewidth',1.5);
hold on 
plot(time(1:end-H),mu_p_hist(:,2),'b-','linewidth',1.5);
legend('True','Estimate')
ylabel('$\Delta_y$','interpreter','latex');
set(gca,'fontsize',16,'Fontname','times new roman');

subplot(3,1,3)
hold on ;grid on;box on
plot(time(1:end-H),target_p_hist(:,3),'r--','linewidth',1.5);
hold on 
plot(time(1:end-H),mu_p_hist(:,3),'b-','linewidth',1.5);
legend('True','Estimate')
xlabel('Times/s');
ylabel('$\Delta_z$','interpreter','latex');
set(gca,'fontsize',16,'Fontname','times new roman');



% std_sd = sqrt(var_sd_hist);      
% 
% figure('Name','GP estimation result-short term ')
% subplot(3,1,1)
% hold on ;grid on;box on
% % patch([time(1:end-H), fliplr(time(1:end-H))],[mu_sp_hist(:,1)-2*std_sd(:,1); flipud(mu_sp_hist(:,1)-2*std_sd(:,1))],[236 241 249]/255, 'EdgeColor', 'none');hold on
% plot(time(1:end-H),target_sp_hist(:,1),'r--','linewidth',1.5);
% hold on 
% plot(time(1:end-H),mu_sp_hist(:,1),'b-','linewidth',1.5);
% legend('True','Estimate')
% set(gca,'fontsize',16,'Fontname','times new roman');
% ylabel('$\Delta_{sx}$','interpreter','latex');
% 
% 
% subplot(3,1,2)
% hold on ;grid on;box on
% plot(time(1:end-H),target_sp_hist(:,2),'r--','linewidth',1.5);
% hold on 
% plot(time(1:end-H),mu_sp_hist(:,2),'b-','linewidth',1.5);
% legend('True','Estimate')
% set(gca,'fontsize',16,'Fontname','times new roman');
% ylabel('$\Delta_{sy}$','interpreter','latex');
% 
% subplot(3,1,3)
% hold on ;grid on;box on
% plot(time(1:end-H),target_sp_hist(:,3),'r--','linewidth',1.5);
% hold on 
% plot(time(1:end-H),mu_sp_hist(:,3),'b-','linewidth',1.5);
% legend('True','Estimate')
% xlabel('Times/s');
% ylabel('Z ');
% set(gca,'fontsize',16,'Fontname','times new roman');
% ylabel('$\Delta_{sz}$','interpreter','latex');

averCost  = sum(cost_trans)/length(cost_trans)

 save('lGPMPC_without_update.mat','x','u','cost_trans','mu_p_hist','MSE','target_p_hist','mu_p_hist');