load para.mat

xd = para.path;
xd =xd';

load baselineMPC.mat
x0 = x;
baseline_MSE = MSE

load lGPMPC_without_update.mat
x1 = x;
lGP_MSE = MSE

% load lGPMPC_with_update.mat
load res_oGP_stack.mat
res = res_oGP_stack{6,2};
x2 = res.x;
oGP_MSE = res.MSE

load res_stack.mat
res = res_stack{4}; %3
DGP_1_MSE = res.MSE

load res_2nd_stack.mat
res = res_2nd_stack{4}; %2
x3  = res.x;

DGP_2_MSE = res.MSE

%% XY plane -comparison
H=5;
green = [167,202,122]/255;
red = [216 56 58]/255;
yellow = [244,205,104]/255;
darkred = [255,122,66]/255;
blue = [47 127 193]/255;


c=para.time;

figure;
ax1 = subplot(1,4,1);
colormap('winter');
cmap = colormap;
c = round(1+(size(cmap,1)-1)*(c - min(c))/(max(c)-min(c)));
plot(x0(1,1:end-H+1),x0(2,1:end-H+1),'color',blue,'linewidth',2)
grid on; hold on; box on 
plot(xd(1,1:end-H+1),xd(2,1:end-H+1),'-.','linewidth',2,'color',red);
plot(x0(1,1),x0(2,1),'d','color',darkred,'Markersize',20,'MarkerFaceColor',darkred);hold on 
plot(xd(1,end-H+1),xd(2,end-H+1),'p','color',yellow,'Markersize',20,'MarkerFaceColor',yellow);hold on 
plot(x0(1,end-H+1),x0(2,end-H+1),'p','color',green,'Markersize',20,'MarkerFaceColor',green);hold on 
for k = 1:(length(x0)-1-H)
    line(x0(1,k:k+1),x0(2,k:k+1),c(k:k+1),'color',cmap(c(k),:),'linewidth',2)
end
% colorbar
axis([-2.5 2.5 -2.5 2.5]);
legend('Actual','Ref','Initial','Final-desired','Final-actual');
xlabel('X/m');
ylabel('Y/m');
set(gca,'fontsize',16,'Fontname','times new roman');


ax2 =subplot(1,4,2);

% colormap('winter');
% cmap = colormap;
% c = round(1+(size(cmap,1)-1)*(c - min(c))/(max(c)-min(c)));
plot(x1(1,1:end-H+1),x1(2,1:end-H+1),'color',blue,'linewidth',2)
grid on; hold on; box on 
plot(xd(1,1:end-H+1),xd(2,1:end-H+1),'-.','linewidth',2,'color',red);
for k = 1:(length(x0)-1-H)
    line(x1(1,k:k+1),x1(2,k:k+1),c(k:k+1),'color',cmap(c(k),:),'linewidth',2)
end
% colorbar
axis([-2.5 2.5 -2.5 2.5]);
plot(x1(1,1),x1(2,1),'d','color',darkred,'Markersize',20,'MarkerFaceColor',darkred);hold on 
plot(xd(1,end-H+1),xd(2,end-H+1),'p','color',yellow,'Markersize',20,'MarkerFaceColor',yellow);hold on 
plot(x1(1,end-H+1),x1(2,end-H+1),'p','color',green,'Markersize',20,'MarkerFaceColor',green);hold on 
set(gca,'fontsize',16,'Fontname','times new roman');
xlabel('X/m');
ylabel('Y/m');


ax3 = subplot(1,4,3);
% colormap('winter');
% cmap = colormap;
% c = round(1+(size(cmap,1)-1)*(c - min(c))/(max(c)-min(c)));
plot(x2(1,1:end-H+1),x2(2,1:end-H+1),'color',blue,'linewidth',2)
grid on; hold on; box on 
plot(xd(1,1:end-H+1),xd(2,1:end-H+1),'-.','linewidth',2,'color',red);
for k = 1:(length(x0)-1-H)
    line(x2(1,k:k+1),x2(2,k:k+1),c(k:k+1),'color',cmap(c(k),:),'linewidth',2)
end
% colorbar
axis([-2.5 2.5 -2.5 2.5]);
plot(x2(1,1),x2(2,1),'d','color',darkred,'Markersize',20,'MarkerFaceColor',darkred);hold on 
plot(xd(1,end-H+1),xd(2,end-H+1),'p','color',yellow,'Markersize',20,'MarkerFaceColor',yellow);hold on 
plot(x2(1,end-H+1),x2(2,end-H+1),'p','color',green,'Markersize',20,'MarkerFaceColor',green);hold on 
set(gca,'fontsize',16,'Fontname','times new roman');
xlabel('X/m');
ylabel('Y/m');


ax4 =subplot(1,4,4);
colormap('winter');
cmap = colormap;
c = round(1+(size(cmap,1)-1)*(c - min(c))/(max(c)-min(c)));
plot(x3(1,1:end-H+1),x3(2,1:end-H+1),'color',blue,'linewidth',2)
grid on; hold on; box on 
plot(xd(1,1:end-H+1),xd(2,1:end-H+1),'-.','linewidth',2,'color',red);
for k = 1:(length(x0)-1-H)
    line(x3(1,k:k+1),x3(2,k:k+1),c(k:k+1),'color',cmap(c(k),:),'linewidth',2)
end
colorbar
caxis([0 para.T]);
axis([-2.5 2.5 -2.5 2.5]);
plot(x3(1,1),x3(2,1),'d','color',darkred,'Markersize',20,'MarkerFaceColor',darkred);hold on 
plot(xd(1,end-H+1),xd(2,end-H+1),'p','color',yellow,'Markersize',20,'MarkerFaceColor',yellow);hold on 
plot(x3(1,end-H+1),x3(2,end-H+1),'p','color',green,'Markersize',20,'MarkerFaceColor',green);hold on 
set(gca,'fontsize',16,'Fontname','times new roman');
xlabel('X/m');
ylabel('Y/m');


ax4.Position([3 4]) = ax1.Position([3 4]);



%% DGP estimation 
load res_2nd_stack.mat

res  = res_2nd_stack{4};
time = res.para.time;
time = time';
target_p_hist = res.target_p_hist;
mu_p_hist = res.mu_p_hist;

target_sp_hist = res.target_sp_hist;
mu_sp_hist = res.mu_sp_hist;


figure('Name','GP estimation result-total ')
subplot(3,1,1)
hold on ;grid on;box on
plot(time(1:end-H),target_p_hist(:,1),'--','color',darkred,'linewidth',3);
hold on 
plot(time(1:end-H),mu_p_hist(:,1),'-','color',blue,'linewidth',3);
legend('True','DGP-Estimation')
ylabel('$\Delta_x$','interpreter','latex');
set(gca,'fontsize',16,'Fontname','times new roman');


subplot(3,1,2)
hold on ;grid on;box on
plot(time(1:end-H),target_p_hist(:,2),'--','color',darkred,'linewidth',3);
hold on 
plot(time(1:end-H),mu_p_hist(:,2),'-','color',blue,'linewidth',3);
% legend('True','Estimate')
ylabel('$\Delta_y$','interpreter','latex');
set(gca,'fontsize',16,'Fontname','times new roman');

subplot(3,1,3)
hold on ;grid on;box on
plot(time(1:end-H),target_p_hist(:,3),'--','color',darkred,'linewidth',3);
hold on 
plot(time(1:end-H),mu_p_hist(:,3),'-','color',blue,'linewidth',3);
% legend('True','Estimate')
xlabel('Time/s');
ylabel('$\Delta_z$','interpreter','latex');
set(gca,'fontsize',16,'Fontname','times new roman');

% figure('Name','GP estimation result-short term ')
% subplot(3,1,1)
% hold on ;grid on;box on
% plot(time(1:end-H),target_sp_hist(:,1),'r--','linewidth',1.5);
% hold on 
% plot(time(1:end-H),mu_sp_hist(:,1),'b-','linewidth',1.5);
% legend('True','Estimate')
% ylabel('$\Delta_x$','interpreter','latex');
% set(gca,'fontsize',16,'Fontname','times new roman');
% 
% 
% subplot(3,1,2)
% hold on ;grid on;box on
% plot(time(1:end-H),target_sp_hist(:,2),'r--','linewidth',1.5);
% hold on 
% plot(time(1:end-H),mu_sp_hist(:,2),'b-','linewidth',1.5);
% legend('True','Estimate')
% ylabel('$\Delta_y$','interpreter','latex');
% set(gca,'fontsize',16,'Fontname','times new roman');
% 
% subplot(3,1,3)
% hold on ;grid on;box on
% plot(time(1:end-H),target_sp_hist(:,3),'r--','linewidth',1.5);
% hold on 
% plot(time(1:end-H),mu_sp_hist(:,3),'b-','linewidth',1.5);
% legend('True','Estimate')
% xlabel('Times/s');
% ylabel('$\Delta_z$','interpreter','latex');
% set(gca,'fontsize',16,'Fontname','times new roman');


%%
res = res_2nd_stack{2};
x  = res.x;
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

%% 