clc;clear;

Max_iter = 25;

res_stack = cell(Max_iter,1);
idx = [];

for i = 1: Max_iter 

   MSE =  DGP_offline_SVGP_MC;

    try
       res = DGPMPC_for_MC;
       
    catch WARN
        warning('Not converged!');
        idx = [idx,i];
        continue;
     
    end

    res_stack{i,1}  = res;

end


%%
x = res_stack{5}.x;
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
