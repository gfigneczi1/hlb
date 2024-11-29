clc;clear;

load res_varGP_1.mat
load res_varGP_2.mat
load res_DGP_1.mat
load res_DGP_2.mat
load dic_GP.mat

fillColor = [236 241 249]/255;

err_varGP_1 = abs(ff'-mustar_varGP_1);
err_varGP_2 = abs(ff'-mustar_varGP_2);
err_DGP_1 = abs(ff'-mustar_DGP_1);
err_DGP_2= abs(ff'-mustar_DGP_2);
err_dicGP_1 = abs(ff'-mustar);
err_dicGP_2 = abs(ff'-mustar_2nd);


blue = [21, 123, 190]/255;
orange = [190, 88, 21]/255;
purple = [128,60,129]/255;
% green = [39, 190, 21]/255;

green = [167,202,122]/255;
red = [216 56 58]/255;
yellow = [244,205,104]/255;
darkred = [255,122,66]/255;
darkblue = [47 127 193]/255;



%%
figure

subplot(3,2,1)
hold on; grid on;box on;
plot(Xtest, ff', '-.','color',darkred, 'LineWidth', 2) % 
plot(Xtest, mustar, '-', 'color',darkblue, 'LineWidth', 2) % mean predictions in blue
legend('True','GP-Estimation');
xlabel('SOGP, 1st iter.');
ylabel('$y$','interpreter','latex');
set(gca, 'fontsize', 16,'fontname','times new roman');

subplot(3,2,2)
hold on; grid on;box on;
plot(Xtest, ff', '-.','color',darkred, 'LineWidth', 2) % 
plot(Xtest, mustar_2nd, '-', 'color',darkblue, 'LineWidth', 2) % mean predictions in blue
xlabel('SOGP, 2nd iter.');
set(gca, 'fontsize', 16,'fontname','times new roman');

subplot(3,2,3)
hold on; grid on;box on;
plot(Xtest, ff', '-.','color',darkred, 'LineWidth', 2) % 
plot(Xtest, mustar_varGP_1, '-', 'color',darkblue,'LineWidth', 2) % mean predictions in blue
xlabel('ROGP, 1st iter.');
ylabel('$y$','interpreter','latex');
set(gca, 'fontsize', 16,'fontname','times new roman');

subplot(3,2,4)
hold on; grid on;box on;
plot(Xtest, ff', '-.','color',darkred, 'LineWidth', 2) % 
plot(Xtest, mustar_varGP_2, '-', 'color',darkblue, 'LineWidth', 2) % mean predictions in blue
xlabel('ROGP, 2nd iter.');
set(gca, 'fontsize', 16,'fontname','times new roman');



subplot(3,2,5)
hold on; grid on;box on;
plot(Xtest, ff', '-.','color',darkred, 'LineWidth', 2) % 
plot(Xtest, mustar_DGP_1, '-', 'color',darkblue, 'LineWidth', 2) % mean predictions in blue
xlabel('DGP, 1st iter.');
ylabel('$y$','interpreter','latex');
set(gca, 'fontsize', 16,'fontname','times new roman');

subplot(3,2,6)

hold on; grid on;box on;
plot(Xtest, ff', '-.','color',darkred, 'LineWidth', 2) % 
plot(Xtest, mustar_DGP_2, '-', 'color',darkblue, 'LineWidth', 2) % mean predictions in blue
xlabel('DGP, 2nd iter.');
set(gca, 'fontsize', 16,'fontname','times new roman');

%%


blue = [21, 123, 190]/255;
orange = [190, 88, 21]/255;
purple = [128,60,129]/255;
green = [39, 190, 21]/255;


figure('Name','GP estimation error','Units','centimeter','Position',[5,5,15,15]);
subplot(2,1,1)
hold on ;grid off;box on    
f11 = plot_shaded(Xtest,err_varGP_1,'color',blue); 
f12 = plot_shaded(Xtest,err_DGP_1,'color',orange); 
f13 = plot_shaded(Xtest,err_dicGP_1,'color',purple);
legend([f11 f12 f13],'1','2','3');
hx = ylabel('GPR Error-1st Iter.');
set(gca,'fontsize',18,'Fontname','times new roman');

subplot(2,1,2)
hold on ;grid off;box on    
f21 = plot_shaded(Xtest,err_varGP_2,'color',blue); 
f22 = plot_shaded(Xtest,err_DGP_2,'color',orange); 
f23 = plot_shaded(Xtest,err_dicGP_2,'color',purple);
legend([f21 f22 f23],'1','2','3');
hx = ylabel('GPR Error-2nd Iter.');
set(gca,'fontsize',18,'Fontname','times new roman');

