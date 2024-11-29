function [x_ref,y_ref,z_ref,vx_ref,vy_ref,vz_ref,ax_ref,ay_ref,az_ref,psi_ref,refTrajectory] = ref_fix_point(time)

x_ref  = 1*ones(1,length(time));
vx_ref = 0*ones(1,length(time));
ax_ref = 0*ones(1,length(time));

y_ref  = 2*ones(1,length(time));
vy_ref = 0*ones(1,length(time));
ay_ref = 0*ones(1,length(time));


z_ref  = 3*ones(1,length(time));
vz_ref = 0*ones(1,length(time));
az_ref = 0*ones(1,length(time));

psi_ref = atan2(vy_ref, vx_ref);

refTrajectory = [x_ref' y_ref' z_ref' vx_ref' vy_ref' vz_ref' ax_ref' ay_ref' az_ref' psi_ref'];