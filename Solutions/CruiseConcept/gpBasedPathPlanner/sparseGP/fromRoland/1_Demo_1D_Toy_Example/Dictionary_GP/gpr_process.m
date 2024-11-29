function [model] = gpr_process(data,y)

 global noise sigma

 %create initial GP model 
    BV = data;
    obs = y';
    current_size = size(data,2);
    
    noise_x = noise + 1; % compute noise param
    Q = y/noise_x;
    C = -1/noise_x;    
    K = kernel_o(data,data,sigma);  % only 1 dimension
    K = K + noise*eye(current_size);
    alpha = y/noise_x;

    model.BV = BV;
    model.obs = obs;
    model.current_size= current_size;
    model.Q = Q;
    model.C = C;
    model.K = K;
    model.alpha = alpha;

end

