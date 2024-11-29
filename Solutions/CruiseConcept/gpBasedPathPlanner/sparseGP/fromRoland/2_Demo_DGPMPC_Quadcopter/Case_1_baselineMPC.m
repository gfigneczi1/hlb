%% Data collection phase
% constant wind added in dynamics
% xdot = Ax+Bu , nominal model
% collect data using the random trajectory

%% Written by Yuhan Liu @Netherlands, July 2020.
clc;clear all; close all
global winds flag_wind

warning off
addpath('util');
addpath('odefiles');
addpath('GP_related');
addpath('MPCsolver');
addpath('varGP');
addpath('pilco-matlab-master');

set_parameters;
init_LTI;

flag_wind = 1; % if there is wind disturbances or not
winds = generate_wind_rough(para,3)*flag_wind;

figure;
plot(winds,'linewidth',2);
%% Dyn settings

para.trans.A  = LTI_trans.A;
para.trans.B  = LTI_trans.B;  % discete 

para.rot.A  =  LTI_rot.A;
para.rot.B  =  LTI_rot.B;  % discete 
% para.gpmodel = gpmodel;


[~,K_rot,~] = idare(para.rot.A,para.rot.B,para.rot.Q,para.rot.R);
kp = diag([2,2,2]);
kd = diag([0.5,0.5,0.5]);

% attitude control
kp_roll_pitch_yaw = [6.0, 6.0,4.0];
% kd_roll_pitch_yaw = [20.0 20.0 20.0]; 
kd_roll_pitch_yaw = [8.0 8.0 8.0];
%%
ref = generate_Ref(para,2);  % 1-->fixed point, 2-->helix, 3-->'8' figure, 4-->random
xd = ref.x(1:6,:);
para.path  = xd';
figure;
plot3(xd(1,:),xd(2,:),xd(3,:));

xtrans0 = xd(:,1);
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

    tic;
     for i = 1: para.N-para.trans.H

        %-----------------------------------------------------------%
        %   Compute the optimal control with in the horizon         %
        %-----------------------------------------------------------%
          
           
          [u_opt_trans,cost_trans(i)]  = mpcCtrl_trans_noComp(xtrans(:,i),time(i),para,i); % Basic MPC 
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
        
        if flag_wind == 1
            [drag_force,drag_torque] = get_delta(i,x(:,i),u(:,i),LTI,para,i);
             drag_force_hist(i,:)    = drag_force';
             drag_torque_hist(i,:)   = drag_torque';
        end
        
        xnew = ode4(@dynsys_Plin_Anonlin_wind,[i-1,i]*dT,x(:,i),u(:,i),LTI,para,i);
        x(:,i+1) = xnew(end,:)';
        xtrans(:,i+1) =  x(1:6,i+1);
        xrot(:,i+1)   =  x(7:12,i+1);
    


     end
        toc

%%
H = para.trans.H;


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

% figure;
% plot(para.time(1,1:end-H+1),x(9,1:end-H+1),'b','linewidth',1.5);
% hold on
% plot(para.time(1,1:end-H+1),ref.psi(1:end-H+1),'r--','linewidth',1.5);
% legend('psi','psi-ref');
% xlabel('Time/s');
% 
% figure;
% subplot(2,1,1)
% plot(para.time(1,1:end-H+1),x(7,1:end-H+1),'b','linewidth',1.5); hold on
% legend('phi','phi-ref');
% 
% subplot(2,1,2)
% plot(para.time(1,1:end-H+1),x(8,1:end-H+1),'b','linewidth',1.5); hold on
% legend('theta','theta-ref');
% xlabel('Time/s');

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



green = [167,202,122]/255;
red = [216 56 58]/255;
yellow = [244,205,104]/255;
darkred = [255,122,66]/255;
blue = [47 127 193]/255;
% [95 151 210]/255;


wind_direction =winds(1,:)/norm(winds(1,:))/2;
figure;
plot3(x(1,1:end-H+1),x(2,1:end-H+1),x(3,1:end-H+1),'color',blue,'linewidth',2)
grid on; hold on; box on 
plot3(xd(1,1:end-H+1),xd(2,1:end-H+1),xd(3,1:end-H+1),'-.','linewidth',2,'color',red);
plot3(x(1,1),x(2,1),x(3,1),'o','color',darkred,'Markersize',20,'MarkerFaceColor',darkred);hold on 
plot3(xd(1,end-H+1),xd(2,end-H+1),xd(3,end-H+1),'p','color',yellow,'Markersize',20,'MarkerFaceColor',yellow);hold on 
plot3(x(1,end-H+1),x(2,end-H+1),x(3,end-H+1),'p','color',green,'Markersize',20,'MarkerFaceColor',green);hold on 
% quiver3(0.2,0.5,3,wind_direction(1),wind_direction(2),wind_direction(3),'linewidth',2,'color','m','MaxHeadSize',0.26);
legend('Actual','Ref','Initial','Final-desired','Final-actual');
xlabel('X,m');
ylabel('Y,m');
zlabel('Z,m');
view(60,10);
set(gca,'fontsize',16,'Fontname','times new roman');


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


% if flag_wind == 1
%     
%     figure('Name','Drag force');
%     plot(time(1:end-H),drag_force_hist(:,1),'b-','linewidth',1.5);hold on ;box on; grid on
%     plot(time(1:end-H),drag_force_hist(:,2),'r--','linewidth',1.5);
%     plot(time(1:end-H),drag_force_hist(:,3),'k:','linewidth',1.5);
%     legend('x','y','z');
%     xlabel('Time/s');
%     title('Drag force');
%     set(gca,'Fontsize',16);
% 
%     figure('Name','Drag torque');
%     plot(time(1:end-H),drag_torque_hist(:,1),'b-','linewidth',1.5);hold on ;box on; grid on
%     plot(time(1:end-H),drag_torque_hist(:,2),'r--','linewidth',1.5);
%     plot(time(1:end-H),drag_torque_hist(:,3),'k:','linewidth',1.5);
%     legend('x','y','z');
%     title('Drag torque');
%     xlabel('Time/s');
%     set(gca,'Fontsize',16);
% 
% end

Err_x = sum((x(1,1:end-H)-xd(1,1:end-H)).^2);
MSE_x = Err_x/(N-H);

Err_y = sum((x(2,1:end-H)-xd(2,1:end-H)).^2);
MSE_y = Err_y/(N-H);

Err_z = sum((x(3,1:end-H)-xd(3,1:end-H)).^2);
MSE_z = Err_z/(N-H);

MSE = [MSE_x MSE_y MSE_z]
  
save('baselineMPC.mat','x','u','cost_trans','MSE','state_dot_hist','LTI');


if flag_wind == 1
    
    figure('Name','Drag force');
    plot(time(1:end-H),drag_force_hist(:,1),'b-','linewidth',1.5);hold on ;box on; grid on
    plot(time(1:end-H),drag_force_hist(:,2),'r--','linewidth',1.5);
    plot(time(1:end-H),drag_force_hist(:,3),'k:','linewidth',1.5);
    legend('x','y','z');
    xlabel('Time/s');
    title('Drag force');
    set(gca,'Fontsize',16);

    figure('Name','Drag torque');
    plot(time(1:end-H),drag_torque_hist(:,1),'b-','linewidth',1.5);hold on ;box on; grid on
    plot(time(1:end-H),drag_torque_hist(:,2),'r--','linewidth',1.5);
    plot(time(1:end-H),drag_torque_hist(:,3),'k:','linewidth',1.5);
    legend('x','y','z');
    title('Drag torque');
    xlabel('Time/s');
    set(gca,'Fontsize',16);

end

        
        averCost  = sum(cost_trans)/length(cost_trans)
