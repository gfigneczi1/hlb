function [mu_xk,var_xk] = multi_pred_ahead(mu_x0,var_x0,uk,t,mpcobj,flag,is_var)   
%------------------------------------------------------------------
        % Propagate mean and covariance multi-step ahead
 %------------------------------------------------------------------
 
   mu_xk  = zeros(mpcobj.n,mpcobj.H+1);
   var_xk = zeros(mpcobj.n,mpcobj.n,mpcobj.H+1);
            
   mu_xk(:,1)    = mu_x0;
   var_xk(:,:,1) = var_x0;
            
        for i = 1:mpcobj.H      % [x1,...,xN]

            [mu_xk(:,i+1),var_xk(:,:,i+1)] = state_pred_xkp1(mu_xk(:,i),var_xk(:,:,i),uk(:,i),t,mpcobj,flag,is_var);  

        end

end






