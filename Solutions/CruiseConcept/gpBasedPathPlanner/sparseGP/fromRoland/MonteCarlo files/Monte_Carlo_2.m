clc;clear;

Max_iter = 25;

res_oGP_stack = cell(Max_iter,2);
idx = [];


for i = 1: Max_iter 

     MSE =  DGP_offline_SVGP_MC;

     try
     res1  = oGPMPC_for_MC();
     res2 = oGPMPC_2nd_for_MC();

     catch WARN
        warning('Not converged!');
        idx = [idx,i];
        continue;
     
    end

    res_oGP_stack{i,1}  = res1;
    res_oGP_stack{i,2}  = res2;





end