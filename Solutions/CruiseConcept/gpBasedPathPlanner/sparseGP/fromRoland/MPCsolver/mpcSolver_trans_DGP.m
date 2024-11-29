function [u_i,f_value] = mpcSolver_trans_DGP(x0,t0)

global mpcobj


persistent u0

if isempty(u0)
    u0 = zeros(mpcobj.m,mpcobj.H);
end
    

% define optimizer settings
options = optimoptions('fmincon',...
                       'Display','off',...
                       'Algorithm', 'interior-point');  % 'sqp','interior-point','active-set'

lbc = zeros(mpcobj.m,1)+mpcobj.lb;
lbc = repmat(lbc,mpcobj.H,1);

ubc = zeros(mpcobj.m,1)+mpcobj.ub;
ubc = repmat(ubc,mpcobj.H,1);


% [u_opt,f_value] = fmincon(@(u)prob_cost(x0,u,t0),u0,[],[],[],[],lbc,ubc,@(u)nonlcon(x0,u,t0),options);
 [u_opt,f_value] = fmincon(@(u)prob_cost(x0,u,t0),u0,[],[],[],[],lbc,ubc,[],options);
 
u0 = u_opt;

u_i = u_opt(:,1);

message = sprintf('NMPC: %%g/%g iterations, Cost = %%20.20g.\n',mpcobj.N);

i_step = t0/mpcobj.dT;
fprintf(message,i_step,f_value);   

end



function [mu_xkp1,var_xkp1] = state_pred(mu_xk, var_xk, uk,t)
    %------------------------------------------------------------------
    %   State prediction (motion model) using exact moment matching
    %
    %       xk+1 = Axk+Buk + delta(xk,uk),  
    %
    %------------------------------------------------------------------
        % calculate discrete time dynamics
global mpcobj
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


      %  DGP prediction
          [mu_dx,var_dx] = DGP_varsgpPredict_x_c_mex(model_x,Smodel_x,mu_xu,var_xu);
          [mu_dy,var_dy] = DGP_varsgpPredict_y_c_mex(model_y,Smodel_y,mu_xu,var_xu);
          [mu_dz,var_dz] = DGP_varsgpPredict_z_c_mex(model_z,Smodel_z,mu_xu,var_xu);

          mu_d  = [mu_dx mu_dy mu_dz]';         
          var_d = diag([var_dx var_dy var_dz]);     % EMM   mu_d  3x1; var_d 3x3


      % predict mean and variance for the state 
     
          mu_xkp1  = f1 + dT*[0;0;0;mu_d] ;    
          var_xkp1 = A*var_xk*A' + dT^2*blkdiag(zeros(3),var_d);
      

end

function [mu_xk,var_xk] = multi_Pred(mu_x0,var_x0,uk,t)   
%------------------------------------------------------------------
        % Propagate mean and covariance multi-step ahead
 %------------------------------------------------------------------
      global mpcobj
   mu_xk  = zeros(mpcobj.n,mpcobj.H+1);
   var_xk = zeros(mpcobj.n,mpcobj.n,mpcobj.H+1);
            
   mu_xk(:,1)    = mu_x0;
   var_xk(:,:,1) = var_x0;
            
        for i = 1:mpcobj.H      % [x1,...,xN]

            [mu_xk(:,i+1),var_xk(:,:,i+1)] = state_pred(mu_xk(:,i),var_xk(:,:,i),uk(:,i),t);  

        end

end



function cost = prob_cost(mu_x0,u,t)   % mu_x0 = x0;  u is the variable waiting for optimization, 1xH
global mpcobj
% initialize the first covariance
var_x0 = zeros(mpcobj.n);
            
% calculate state sequence for given control input sequence and x0
[mu_xk,var_xk] = multi_Pred(mu_x0, var_x0, u,t);    

cost = 0;
path = mpcobj.path;  % n x 
Q = mpcobj.Q;
R = mpcobj.R;
P = mpcobj.P;
H = mpcobj.H;
dT = mpcobj.dT;

 nn = round((t+dT:dT:t+(H+1)*dT)/dT);
 ref = path(nn,:); % H+1 x n
 ref = ref';   % n x H+1 , each column is a desired state for one step

% stage cost accumulation
for i = 1:mpcobj.H
    
    cost = cost + (mu_xk(:,i)-ref(:,i))'*Q*(mu_xk(:,i)-ref(:,i))+ u(:,i)'*R*u(:,i)...
                 +trace(Q*var_xk(:,:,i));
                
    
end

% terminal cost
cost = cost + (mu_xk(:,end)-ref(:,end))'*P*(mu_xk(:,end)-ref(:,end))...
            +trace(P*var_xk(:,:,end));
            

end

function [cineq,ceq] = nonlcon(mu_x0,u,t)

global mpcobj
               
 % init outputs
    ceq = []; 
    cineq = [];

            H = mpcobj.H;
            ng  = 6;
       
         % initialize the first covariance
            var_x0 = zeros(mpcobj.n);
                  
        % init vectors to speedup calculations
%             ceq_h   = zeros(obj.nh, H);
            cineq_gl = zeros(ng, H);
            cineq_gu = zeros(ng, H);
                   
            % calculate state sequence for given control input sequence and x0
            [mu_xk,var_xk] = multi_Pred(mu_x0,var_x0,u,t); 
     
            for i = 1:H
              
                cineq_gl(:,i) = low_con(mu_xk(:,i),var_xk(:,:,i),u,t);
                cineq_gu(:,i) = up_con(mu_xk(:,i),var_xk(:,:,i),u,t);

            end

%             ceq   = ceq_h(:);
            cineq = [cineq_gl(:);cineq_gu(:)];   % column output



end



function cineq_l = low_con(mu,var,u,t)
% x -  6x1 ,var - 6x6 
global mpcobj

con   = mpcobj.con;
x_min = con(1);
y_min = con(2);
z_min = con(2);
vx_min = con(4);
vy_min = con(5);
vz_min = con(6);

sigma = sqrt(diag(var));

c1 = x_min - mu(1) + 2*sigma(1);
c2 = y_min - mu(2) + 2*sigma(2);
c3 = z_min - mu(3) + 2*sigma(3);

c4 = vx_min - mu(4) + 2*sigma(4);
c5 = vy_min - mu(5) + 2*sigma(5);
c6 = vz_min - mu(6) + 2*sigma(6);

cineq_l = [c1 c2 c3 c4 c5 c6]';

end

function cineq_u = up_con(mu,var,u,t)
% x -  6x1 ,var - 6x6 
global mpcobj

con   = mpcobj.con;
x_max = con(7);
y_max = con(8);
z_max = con(9);
vx_max = con(10);
vy_max = con(11);
vz_max = con(12);

sigma = sqrt(diag(var));

c1 = mu(1) + 2*sigma(1) - x_max;
c2 = mu(2) + 2*sigma(2) - y_max;
c3 = mu(3) + 2*sigma(3) - z_max;

c4 = mu(4) + 2*sigma(4) - vx_max;
c5 = mu(5) + 2*sigma(5) - vy_max;
c6 = mu(6) + 2*sigma(6) - vz_max;

cineq_u = [c1 c2 c3 c4 c5 c6]';

end