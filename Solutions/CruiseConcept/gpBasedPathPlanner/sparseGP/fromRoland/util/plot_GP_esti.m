clc;clear;

load lGPMPC_without_update.mat
err_no_update = abs(target_p_hist-mu_p_hist);

load lGPMPC_with_update.mat
load res_oGP_stack.mat
res = res_oGP_stack{2,2};
err_with_update = abs(res.target_p_hist-res.mu_p_hist);

load res_stack.mat
res = res_stack{2};
err_DGP_1st = abs(res.target_p_hist-res.mu_p_hist);

load res_2nd_stack.mat
res = res_2nd_stack{4}; %2
err_DGP_2nd = abs(res.target_p_hist-res.mu_p_hist);

time = res.para.time;

% load GP_esti_no_update;
% err_no_update = abs(target_p_hist-mu_p_hist);
% 
% load GP_esti_with_update;
% err_with_update = abs(target_p_hist-mu_p_hist);
% 
% load DGP_3rd;
% err_DGP = abs(target_p_hist-mu_p_hist);


orange = [0.85 0.33 0.10]; % orange
blue = [0 0.45 0.74];    % blue
purple = [0.49 0.18 0.56]; % purple
yellow =[0.93 0.69 0.13]; % yellow
 % [0.85 0.33 0.10] % orange
% [0 0.45 0.74]    % blue
% [0.49 0.18 0.56] % purple
% [0.93 0.69 0.13] % yellow

% load GP_esti_with_update
H = 5;
figure('Name','GP estimation error ')
subplot(3,1,1)
hold on ;grid on;box on
f11 = plot_shaded(time(1:end-H),err_no_update(:,1),'color',blue); 
hold on 
f12 = plot_shaded(time(1:end-H),err_with_update(:,1),'color',yellow);
% f13 = plot_shaded(time(1:end-H),err_DGP_1st(:,1),'color',purple);
f14 = plot_shaded(time(1:end-H),err_DGP_2nd(:,1),'color',purple);
% legend(f,'Baseline GPMPC','Online GPMPC','DGPMPC');
legend([f11 f12 f14],'LGP','OGP','DGP');
ylabel('GPR Error-x');
set(gca,'fontsize',16,'Fontname','times new roman');



subplot(3,1,2)
hold on ;grid on;box on
plot_shaded(time(1:end-H),err_no_update(:,2),'color',blue); 
hold on 
plot_shaded(time(1:end-H),err_with_update(:,2),'color',yellow);
plot_shaded(time(1:end-H),err_DGP_2nd(:,2),'color',purple);
ylabel('GPR Error-y');
set(gca,'fontsize',16,'Fontname','times new roman');

subplot(3,1,3)
hold on ;grid on;box on
plot_shaded(time(1:end-H),err_no_update(:,3),'color',blue); 
hold on 
plot_shaded(time(1:end-H),err_with_update(:,3),'color',yellow);
plot_shaded(time(1:end-H),err_DGP_2nd(:,3),'color',purple);
yticks([0,0.1]);
xlabel('Time/s');
ylabel('GPR Error-z');
set(gca,'fontsize',16,'Fontname','times new roman');

% figure;
% subplot(3,1,1)
% plot_shaded(time(1:end-H),err_no_update(:,3),'color',[0.93 0.69 0.13]); 
% [0.85 0.33 0.10] % orange
% [0 0.45 0.74]    % blue
% [0.49 0.18 0.56] % purple
% [0.93 0.69 0.13] % yellow

%% plot bar

N1 =  198;
N2 = 396;

for i = 1:N1
   
     err_DGP_01(:,i) = err_DGP_2nd(i,:)';
     err_OGP_01(:,i) = err_with_update(i,:)';

end

MSE_DGP_01_X = sum(err_DGP_01(1,:).^2)/N1;
MSE_DGP_01_Y = sum(err_DGP_01(2,:).^2)/N1;
MSE_DGP_01_Z = sum(err_DGP_01(3,:).^2)/N1;
MSE_DGP_01 = norm([MSE_DGP_01_X MSE_DGP_01_Y MSE_DGP_01_Z])


MSE_OGP_01_X = sum(err_OGP_01(1,:).^2)/N1;
MSE_OGP_01_Y = sum(err_OGP_01(2,:).^2)/N1;
MSE_OGP_01_Z = sum(err_OGP_01(3,:).^2)/N1;
MSE_OGP_01 = norm([MSE_OGP_01_X MSE_OGP_01_Y MSE_OGP_01_Z])



for i = N1+1:N2
    
      err_DGP_12(:,i) = err_DGP_2nd(i,:)';
      err_OGP_12(:,i) = err_with_update(i,:)';
end


MSE_DGP_12_X = sum(err_DGP_12(1,:).^2)/N1;
MSE_DGP_12_Y = sum(err_DGP_12(2,:).^2)/N1;
MSE_DGP_12_Z = sum(err_DGP_12(3,:).^2)/N1;
MSE_DGP_12 = norm([MSE_DGP_12_X MSE_DGP_12_Y MSE_DGP_12_Z])

MSE_OGP_12_X = sum(err_OGP_12(1,:).^2)/N1;
MSE_OGP_12_Y = sum(err_OGP_12(2,:).^2)/N1;
MSE_OGP_12_Z = sum(err_OGP_12(3,:).^2)/N1;
MSE_OGP_12 = norm([MSE_OGP_12_X MSE_OGP_12_Y MSE_OGP_12_Z])


MSE_DGP = [MSE_DGP_12;MSE_DGP_01];
MSE_OGP = [MSE_OGP_12;MSE_OGP_01];


MSE_res= [MSE_DGP MSE_OGP];
width = 0.8;
figure;hold on;box on;grid on
h = bar3(MSE_res,width);
zlabel('MSE of GP Estimation');
set(gca,'yticklabel',{'','10-20s','','0-10s',''});
set(gca,'xticklabel',{'','DGP','','OGP','',''});
set(h(2),'facecolor',blue);
set(h(1),'facecolor',yellow);
set(gca,'fontsize',16,'fontname','times new roman');


