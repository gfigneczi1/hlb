clc;clear;

maxIter = 10;


for i = 1: maxIter
   tol = 0.1*rand; band = 1*rand;
    time_0(i,:) = MC_varGP();
    time_1(i,:) = MC_DGP();
     time_2(i,:) = MC_SOGP(tol,band);
end





%% plot result


s_blue=[151, 191, 224]/255;
d_blue =  [47 127 193]/255;
d_red = [255,122,66]/255;
s_red = [255, 191, 161]/255;
green = [167,202,122]/255;
purple = [128,60,129]/255;
s_purple = [222, 185, 223]/255;
s_green = [207, 226, 182]/255;


aver_time0 = mean(time_0);
aver_time1 = mean(time_1);
aver_time2 = mean(time_2);

reduction_varGP = (aver_time2 - aver_time0)./aver_time2*100;
reduction_DGP = (aver_time2 - aver_time1)./aver_time2*100;

figure; 
subplot(2,1,1)10
hold on; box on;grid on 
for i = 1: maxIter
plot(time_0(i,:),'color',s_blue,'linewidth',0.5);
plot(time_1(i,:),'color',s_red,'linewidth',0.5);
 plot(time_2(i,:),'color',s_purple ,'linewidth',0.5);
% plot(time_2(i,:),'color',s_green ,'linewidth',0.5);


end
f1=plot(aver_time0,'color',d_blue,'linewidth',3);
f2 = plot(aver_time1,'color',d_red,'linewidth',3);
f3 = plot(aver_time2,'color',purple,'linewidth',3);
ylabel('Comp. time (s)');
legend([f1 f2 f3],'ROGP','DGP','SOGP');
axis([0 305 0 0.01]);
set(gca,'fontsize',18,'Fontname','times new roman');

% 
% figure;
% hold on 
% for i = 1: maxIter
% plot(time_2(i,:),'color',green);
% end

for i = 1:length(reduction_varGP)
    if reduction_varGP(i)<0
        reduction_varGP(i) = 0;
    end

    if reduction_DGP(i)<0
        reduction_DGP(i) = 0;
    end
end

subplot(2,1,2)
plot(reduction_varGP,'color',d_blue,'linewidth',2); hold on; grid on 
plot(reduction_DGP,'color',d_red,'linewidth',2);
ylabel('Comp. reduction (%)');
xlabel('Step');

axis([0 305 0 100]);
set(gca,'fontsize',18,'Fontname','times new roman');
