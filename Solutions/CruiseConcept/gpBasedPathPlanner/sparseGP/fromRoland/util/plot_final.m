  %% plot the result   

%     global winds 
para.T  = 20;       % simulation time (seconds)
% para.dT = 0.05;      % step size (seconds)
 para.dT = 0.1;
time = 0:para.dT:para.T;  % time span
[x_ref,y_ref,z_ref,vx_ref,vy_ref,vz_ref,ax_ref,ay_ref,az_ref,psi_ref,~] = ref_helix(time,2,3,4);
%     wind_direction =winds(1,:)/norm(winds(1,:))/2;
    x_flight = x(1,:);
    y_flight = x(2,:);
    z_flight = x(3,:);
    
    figure;
    plot3(x_ref, y_ref, z_ref, 'r--'); 
    grid on; hold on;
    plot3(x_flight, y_flight, z_flight, 'b','linewidth',1.5);hold on;
    plot3(x_flight(1),y_flight(1),z_flight(1),'ro','Markersize',20,'MarkerFaceColor',[.49 1 .63]);hold on 
    plot3(x_ref(end),y_ref(end),z_ref(end),'p','Markersize',20,'MarkerFaceColor',[.49 1 .63]);hold on 
    plot3(state_hist(end,1),state_hist(end,2),state_hist(end,3),'p','Markersize',20,'MarkerFaceColor',[1 0 0]);hold on 
%     quiver3(0.2,0.5,3,wind_direction(1),wind_direction(2),wind_direction(3),'linewidth',2,'color','m','MaxHeadSize',0.26);
    legend('Ref','Actual','Initial','Final-desired','Final-actual');
%     legend('Ref','Actual');
    xlabel('X,m');
    ylabel('Y,m');
    zlabel('Z,m');
%     title('Flight Trajectory During Data Collection Phase')
    view(60,10);
   set(gca,'fontsize',16,'Fontname','times new roman');
    
    figure;
    subplot(3,2,[1,3,5])
    plot3(x_ref, y_ref, z_ref, 'r--'); grid on; hold on;
    plot3(x_flight(1:end-5), y_flight(1:end-5), z_flight(1:end-5), 'b','linewidth',1.5);hold on;
    plot3(x_flight(1),y_flight(1),z_flight(1),'ro','Markersize',20,'MarkerFaceColor',[.49 1 .63]);hold on 
    plot3(x_ref(end),y_ref(end),z_ref(end),'p','Markersize',20,'MarkerFaceColor',[.49 1 .63]);hold on 
    plot3(x(1,end-5),x(2,end-5),x(3,end-5),'p','Markersize',20,'MarkerFaceColor',[1 0 0]);hold on 
%     quiver3(-0.7,0,1.5,wind_direction(1),wind_direction(2),wind_direction(3),'linewidth',2,'color','m','MaxHeadSize',0.26);
    legend('Ref','Actual','Initial','Final-desired','Final-actual','wind');
    xlabel('X,m');
    ylabel('Y,m');
    zlabel('Z,m');
    view(65,30);
     set(gca,'fontsize',16,'Fontname','times new roman');
    
    subplot(3,2,2)
    plot(time(1:end-5),x_flight(1:end-5),'b-','linewidth',1.5);
    hold on;grid on;
    plot(time(1:end-5),x_ref(1:end-5),'r--','linewidth',1);
    xlabel('Time/ s');
    ylabel('Position-X/ m');
    legend('Actual','Ref');
     set(gca,'fontsize',16,'Fontname','times new roman');
    
    
    subplot(3,2,4)
    plot(time(1:end-5),y_flight(1:end-5),'b-','linewidth',1.5);
    hold on;grid on;
    plot(time(1:end-5),y_ref(1:end-5),'r--','linewidth',1);
    xlabel('Time/s');
    ylabel('Position-Y/ m');
    legend('Actual','Ref');
     set(gca,'fontsize',16,'Fontname','times new roman');
  
    
    subplot(3,2,6)
    plot(time(1:end-5),z_flight(1:end-5),'b-','linewidth',1.5);
    hold on;grid on;
    plot(time(1:end-5),z_ref(1:end-5),'r--','linewidth',1);
    xlabel('Time/s');
    ylabel('Position-Z/ m');
    legend('Actual','Ref');
    set(gca,'fontsize',16,'Fontname','times new roman');
    
    
    
    
%     figure('Name','The reference signal');
%     plot3(x_ref,y_ref,z_ref,'b-','linewidth',1.5)
%     grid on;hold on;box on;
%     plot3(x_ref(1),y_ref(1),z_ref(1),'ro','Markersize',5);
%     xlabel('X,m');
%     ylabel('Y,m');
%     zlabel('Z,m');
%     legend('Ref trajectory')
%     view(65,30);
    


  
    figure('Name','Position X');
    plot(time,x_flight,'b-','linewidth',1.5);
    hold on;grid on;
    plot(time,x_ref,'r--','linewidth',1);
    xlabel('Time, s');
    ylabel('Position-X, m');
    legend('Actual','Ref');
    set(gca,'fontsize',15);
    
    figure('Name','Position Y');
    plot(time,y_flight,'b-','linewidth',1.5);
    hold on;grid on;
    plot(time,y_ref,'r--','linewidth',1);
    xlabel('Time, s');
    ylabel('Position-Y, m');
    legend('Actual','Ref');
    set(gca,'fontsize',15);
    
    figure('Name','Position Z');
    plot(time,z_flight,'b-','linewidth',1.5);
    hold on;grid on;
    plot(time,z_ref,'r--','linewidth',1);
    xlabel('Time, s');
    ylabel('Position-Z, m');
    legend('Actual','Ref');
    set(gca,'fontsize',15);
    
    figure('Name','Position Error');
    err_x = (x_ref - x_flight').^2;
    err_y = (y_ref - y_flight').^2;
    err_z = (z_ref - z_flight').^2;
    
    err_pos = sqrt(err_x + err_y + err_z);
    plot(time, err_pos, 'b','linewidth',1.5); grid on;
    xlabel('Time, s');
    ylabel('Position Error, m');
    
    figure('Name','Position Error of X, Y and Z')
    err_1 = x_ref - x_flight'; 
    err_2 = y_ref - y_flight';
    err_3 = z_ref - z_flight';
    
    plot(time, err_1, 'b','linewidth',1.5); grid on; hold on 
    plot(time, err_2, 'r--','linewidth',1.5);
    plot(time, err_3, 'k:','linewidth',1.5);
    legend('X-axis','Y-axis','Z-axis')
    
    
    
    figure('Name','Yaw-psi');
    quad_psi = state_hist(:,6)';
    plot(time, quad_psi, 'b','linewidth',1.5);
    hold on;  grid on; 
    plot(time, psi_ref, 'r--');
    xlabel('Time, s');
    ylabel('Yaw angle, rad');
    legend('Actual','Ref');
    set(gca,'fontsize',15);
    
    
    figure('Name','Thrust');
    plot(time,F_hist,'b-','linewidth',1.5);
    grid on
    xlabel('Time, s');
    ylabel('Thrust force');
    set(gca,'fontsize',15);

    figure('Name','Torque');
    plot(time,torque_hist(:,1)','b-','linewidth',1.5);
    hold on 
    plot(time,torque_hist(:,2)','r--','linewidth',1.5);
    hold on
    plot(time,torque_hist(:,3)','g.-','linewidth',1.5);
    grid on
    xlabel('Time, s');
    ylabel('Control Torque, N.m');
    legend('tau_x','tau_y','tau_z');
    set(gca,'fontsize',15);
    
%     if if_GP_on == 1
%     figure('Name','GP estimation with variance ')
%     mu1 = mu_p_hist(:,1);
% %     var1 = var_p_hist(:,1); var1 = sqrt(var1);
%     real1 = getDiff_p(state_hist,statedot_hist,R_3_hist,F_hist);
%     f = [mu1+2*var1; flipdim(mu1-2*var1,1)];
%     fill([time'; flipdim(time',1)], f, [7 7 7]/8)
%     hold on;
%     plot(time,mu1,'g')
%     plot(time,real1,'b-')
%     hold on;
%   
%     xlabel('Time')
%     end
%   

norm(err_1)
norm(err_2)
norm(err_3)
%     
%     
% Err_1_pd =  err_1;
% Err_2_pd = err_2;
% Err_3_pd = err_3;
% save('Err_pd.mat','Err_1_pd','Err_2_pd','Err_3_pd');
%     
    
    
    