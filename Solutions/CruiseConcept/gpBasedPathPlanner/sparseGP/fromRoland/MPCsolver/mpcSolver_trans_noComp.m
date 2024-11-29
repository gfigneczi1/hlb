function [u_i,f_value] = mpcSolver_trans_noComp(x0,t0)

global mpcobj
persistent u0

if isempty(u0)
    u0 = zeros(mpcobj.m,mpcobj.H);
end
    

% define optimizer settings
options = optimoptions('fmincon',...
                       'Display','off',...
                       'Algorithm', 'interior-point'); % 'interior-point',... % 'sqp','interior-point','active-set'

lbc = zeros(mpcobj.m,1)-inf;
lbc = repmat(lbc,mpcobj.H,1);

ubc = zeros(mpcobj.m,1)+inf;
ubc = repmat(ubc,mpcobj.H,1);

[u_opt,f_value] = fmincon(@(u)prob_cost(x0,u,t0),u0,[],[],[],[],lbc,ubc,[],options);
 
 u0 = u_opt;

u_i = u_opt(:,1);

message = sprintf('NMPC: %%g/%g iterations, Cost = %%20.20g.\n',mpcobj.N);

fprintf(message,mpcobj.i,f_value);   

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
        A = mpcobj.A;
        B = mpcobj.B;

        dT = mpcobj.dT;
        

        f1 = A*mu_xk + B*uk ;   % part 1


      % combine GP inputs
    
      % evaluate disturbance
        

      % predict mean and variance 
          mu_xkp1  = f1 ;
          var_xkp1 = zeros(6);
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
path = mpcobj.path;  % n x (H+1)
Q = mpcobj.Q;
R = mpcobj.R;
P = mpcobj.P;
H = mpcobj.H;
dT = mpcobj.dT;

 nn = round((t+dT:dT:t+(H+1)*dT)/dT);
 ref = path(nn,:); % (H+1)*2
 ref = ref';

% stage cost accumulation
for i = 1:mpcobj.H
    
    cost = cost + (mu_xk(:,i)-ref(:,i))'*Q*(mu_xk(:,i)-ref(:,i))+ u(:,i)'*R*u(:,i)...
                 +trace(Q*var_xk(:,:,i));
                
    
end

% terminal cost
cost = cost + (mu_xk(:,end)-ref(:,end))'*P*(mu_xk(:,end)-ref(:,end))...
            +trace(P*var_xk(:,:,end));
            

end

function [cineq,ceq] = nonlcon(t0, mu_x0)
global mpcobj
        
ceq = []; 
cineq = [];


end