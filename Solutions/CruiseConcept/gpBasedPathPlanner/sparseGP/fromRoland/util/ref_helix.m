function  [x_ref,y_ref,z_ref,vx_ref,vy_ref,vz_ref,ax_ref,ay_ref,az_ref,psi_ref,refTrajectory]  = ref_helix(time,z1,z2,z3)

omega = 0.75; 
T = time(end);
% 
% n1 = ceil(length(time)*0.3);
% n2 = floor(length(time)*0.3);
% n3 = floor(length(time)*0.4);
d = 2/T;

x_arg   = omega*time;
x_ref   = 2*sin(x_arg);
vx_ref  = 2*omega*cos(x_arg);
ax_ref = -2*omega*omega*sin(x_arg);

y_arg   = omega*time;
y_ref   = 2*cos(y_arg);
vy_ref  = -2*omega*sin(y_arg);
ay_ref = -2*omega*omega*cos(y_arg);


z_ref   = d*time+2;
vz_ref  = d*ones(1,length(time));
az_ref = 0*ones(1,length(time));

% psi_ref = atan2(vy_ref, vx_ref);
 psi_ref = 0*ones(1,length(time));
refTrajectory = [x_ref' y_ref' z_ref' vx_ref' vy_ref' vz_ref' ax_ref' ay_ref' az_ref' psi_ref'];