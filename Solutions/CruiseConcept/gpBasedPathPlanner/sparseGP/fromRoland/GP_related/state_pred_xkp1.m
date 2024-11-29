function [mu_xkp1,var_xkp1] = state_pred_xkp1(mu_xk, var_xk, uk,t,mpcobj,flag,is_var)
    %------------------------------------------------------------------
    %   State prediction (motion model) using exact moment matching
    %
    %       xk+1 = Axk+Buk + delta(xk,uk),  
    %
    %------------------------------------------------------------------
        % calculate discrete time dynamics

        A  = mpcobj.A;
        B  = mpcobj.B;
        dT = mpcobj.dT;

        model_x = mpcobj.model_x;
        model_y = mpcobj.model_y;
        model_z = mpcobj.model_z;
        
        Smodel_x = mpcobj.Smodel_x;
        Smodel_y = mpcobj.Smodel_y;
        Smodel_z = mpcobj.Smodel_z;
        
        xrot     = mpcobj.xrot;

        f1 = A*mu_xk + B*uk ;   % part 1  -----> linear
        

      % combine GP inputs
    
         mu_xu = [xrot(1:3);mu_xk(4:6);uk(3)]';     % test point: [phi theta psi vx vy vz u(3)]' 7x1
        
         var_xu = blkdiag(zeros(3),var_xk(4:6,4:6),0);      % 7x7 blocked covariance matrix 

        % evaluate disturbance

%      % long-term
%           [mul_x, varl_x] = varsgpPredict_post_x(model_x,mu_xu);
%           [mul_y, varl_y] = varsgpPredict_post_y(model_y,mu_xu);
%           [mul_z, varl_z] = varsgpPredict_post_z(model_z,mu_xu);
% 
%           mu_l = [mul_x mul_y mul_z]' ;
%      % short-term
%           [mus_x, vars_x] = varsgpPredict_post_x(Smodel_x,mu_xu);
%           [mus_y, vars_y] = varsgpPredict_post_y(Smodel_y,mu_xu);
%           [mus_z, vars_z] = varsgpPredict_post_z(Smodel_z,mu_xu);
%          
%           mu_s  = [mus_x mus_y mus_z]';
%          
%           mu_d = mu_l + mu_s;

     if flag == 1  % long-term only
          [~,~,mu_dx,var_dx] = varsgpPredict_MM_x(model_x,mu_xu,var_xu);
          [~,~,mu_dy,var_dy] = varsgpPredict_MM_y(model_y,mu_xu,var_xu);
          [~,~,mu_dz,var_dz] = varsgpPredict_MM_z(model_z,mu_xu,var_xu);
     elseif flag == 2 % short term only 
          [~,~,mu_dx,var_dx] = varsgpPredict_MM_x(Smodel_x,mu_xu,var_xu);
          [~,~,mu_dy,var_dy] = varsgpPredict_MM_y(Smodel_y,mu_xu,var_xu);
          [~,~,mu_dz,var_dz] = varsgpPredict_MM_z(Smodel_z,mu_xu,var_xu);
     elseif flag == 3 % dual case

          [mu_dx,var_dx] = DGP_varsgpPredict_x(model_x,Smodel_x,mu_xu,var_xu);
          [mu_dy,var_dy] = DGP_varsgpPredict_y(model_y,Smodel_y,mu_xu,var_xu);
          [mu_dz,var_dz] = DGP_varsgpPredict_z(model_z,Smodel_z,mu_xu,var_xu);
     end


          mu_d  = [mu_dx mu_dy mu_dz]';         
          var_d = diag([var_dx var_dy var_dz]);     % EMM   mu_d  3x1; var_d 3x3


      % predict mean and variance for the state 
     
          mu_xkp1  = f1 + dT*[0;0;0;mu_d] ;    
          var_xkp1 = A*var_xk*A' + dT.^2*blkdiag(zeros(3),var_d);

       if is_var == 1 

          var_xkp1 = A*var_xk*A' + dT.^2*blkdiag(zeros(3),var_d);

       elseif is_var == 0

          var_xkp1 = zeros(6);

       end
end

