%% SVGP

codegen varsgpPredict_MM_x_c -args {model_x,zeros(1,7),zeros(7)}

codegen varsgpPredict_MM_y_c -args {model_y,zeros(1,7),zeros(7)}

codegen varsgpPredict_MM_z_c -args {model_z,zeros(1,7),zeros(7)}

%% DGP
codegen DGP_varsgpPredict_x_2 -args {model_x,Smodel_x,zeros(1,7),zeros(7)}

codegen DGP_varsgpPredict_y_2 -args {model_y,Smodel_y,zeros(1,7),zeros(7)}

codegen DGP_varsgpPredict_z_2 -args {model_z,Smodel_z,zeros(1,7),zeros(7)}

%% GP update

% codegen GP_update_SVGP  -args {zeros(1,7),zeros(1,3),zeros(1,3),Smodel_x,Smodel_y,Smodel_z}


%% test speed

aa = rand(1,7);
bb = rand(7);

s1 = tic;
[mu,s2]=DGP_varsgpPredict_x_c(model_x,Smodel_x,aa,bb)
toc(s1)

ss2 = tic;
[mu,s2]=DGP_varsgpPredict_x_c_mex(model_x,Smodel_x,aa,bb)
toc(ss2)