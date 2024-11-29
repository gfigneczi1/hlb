function  [x_ref,y_ref,z_ref,vx_ref,vy_ref,vz_ref,ax_ref,ay_ref,az_ref,psi_ref,refTrajectory]  = ref_eight_figure(time)

omega_x = 0.8; omega_y = 0.4; omega_z = 0.4;
x_arg   = omega_x*time;
x_ref   = sin(x_arg);
vx_ref  = omega_x*cos(x_arg);
ax_ref = -omega_x*omega_x*sin(x_arg);

y_arg   = omega_y*time;
y_ref   = cos(y_arg);
vy_ref  = -omega_y*sin(y_arg);
ay_ref = -omega_y*omega_y*cos(y_arg);

z_arg   = omega_z*time;
z_ref   = cos(z_arg)+1;
vz_ref  = -omega_z*sin(z_arg);
az_ref = -omega_z*omega_z*cos(z_arg);

psi_ref = atan2(vy_ref, vx_ref);
refTrajectory = [x_ref' y_ref' z_ref' vx_ref' vy_ref' vz_ref' ax_ref' ay_ref' az_ref' psi_ref'];