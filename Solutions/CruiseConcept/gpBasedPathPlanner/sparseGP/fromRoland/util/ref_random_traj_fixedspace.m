function  [x_ref,y_ref,z_ref,vx_ref,vy_ref,vz_ref,ax_ref,ay_ref,az_ref,psi_ref,refTrajectory]  = ref_random_traj_fixedspace(time)

omega_x1 = 1; omega_x2 = 0.33; 
omega_y1 = 1.1; omega_y2 = 0.11;
omega_z1 = 0.33; omega_z2 = 0.76;

ampx1 = 1.5; ampx2 =1.5;
% ampy1 =  0.5; ampy2 = 0.5;
ampy1 =  1; ampy2 =1;
ampz1 = 1; ampz2 = 1;

x_arg1   = omega_x1*time;
x_arg2   = omega_x2*time;
x_ref    = ampx1*sin(x_arg1)  +  ampx2*sin(x_arg2) ;
vx_ref   = ampx1*omega_x1*cos(x_arg1)  +  ampx2*omega_x2*cos(x_arg2);
ax_ref   = -ampx1*omega_x1*omega_x1*sin(x_arg1)  - ampx2*omega_x2*omega_x2*sin(x_arg2);

y_arg1   = omega_y1*time+pi/2;
y_arg2   = omega_y2*time;
% y_ref    = ampy1*sin(y_arg1)  +  ampy2*sin(y_arg2)+1;
y_ref    = ampy1*sin(y_arg1)  +  ampy2*sin(y_arg2);
vy_ref   = ampy1*omega_y1*cos(y_arg1)  + ampy2*omega_y2*cos(y_arg2);
ay_ref   = -ampy1*omega_y1*omega_y1*sin(y_arg1)  -  ampy2*omega_y2*omega_y2*sin(y_arg2);

z_arg1   = omega_z1*time;
z_arg2   = omega_z2*time;
z_ref    = ampz1*sin(z_arg1)  +  ampz2*sin(z_arg2)+3;
vz_ref  = ampz1*omega_z1*cos(z_arg1)  +  ampz2*omega_z2*cos(z_arg2);
az_ref = -ampz1*omega_z1*omega_z1*sin(z_arg1)  - ampz2*omega_z2*omega_z2*sin(z_arg2);

psi_ref = atan2(vy_ref, vx_ref);
refTrajectory = [x_ref' y_ref' z_ref' vx_ref' vy_ref' vz_ref' ax_ref' ay_ref' az_ref' psi_ref'];