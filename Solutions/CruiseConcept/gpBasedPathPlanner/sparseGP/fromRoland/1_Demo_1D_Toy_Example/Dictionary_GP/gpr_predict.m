function [mu,s2] = gpr_predict(model,x)

 global sigma noise

  BV = model.BV;
  alpha = model.alpha;
  C = model.C;

    k = kernel_o(x,BV,sigma)';
    mu = k'*alpha;
    s2 = kernel_o(x,x,sigma) + k'*C*k;  


end

 