function ref = generate_Ref(para,ref_traj)

time = para.time;

switch ref_traj
       case 1
           % fix point
  [x_ref,y_ref,z_ref,vx_ref,vy_ref,vz_ref,ax_ref,ay_ref,az_ref,psi_ref,~] = ref_fix_point(time);
       case 2
           % helix 
  [x_ref,y_ref,z_ref,vx_ref,vy_ref,vz_ref,ax_ref,ay_ref,az_ref,psi_ref,~] = ref_helix(time,2,3,4);
       case 3
           % 'eight'-figure
  [x_ref,y_ref,z_ref,vx_ref,vy_ref,vz_ref,ax_ref,ay_ref,az_ref,psi_ref,~] = ref_eight_figure(time);
       case 4
           % pseudo-random
  [x_ref,y_ref,z_ref,vx_ref,vy_ref,vz_ref,ax_ref,ay_ref,az_ref,psi_ref,~] = ref_random_traj_fixedspace(time);
end

% ref.x = x_ref;
% ref.y = y_ref;
% ref.z = z_ref;
% 
% ref.vx = vx_ref;
% ref.vy = vy_ref;
% ref.vz = vz_ref;

% ref.ax = ax_ref;
% ref.ay = ay_ref;
% ref.az = az_ref;

ref.psi = psi_ref;

ref.x = [x_ref;y_ref;z_ref;vx_ref;vy_ref;vz_ref;ax_ref;ay_ref;az_ref];
% ref.x.rot   = []
