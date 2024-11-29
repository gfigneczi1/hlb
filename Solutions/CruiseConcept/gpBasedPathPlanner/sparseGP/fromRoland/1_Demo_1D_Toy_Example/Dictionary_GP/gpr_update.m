function model_1 = gpr_update(model,x,y,tol,max_points)

global sigma noise

  BV = model.BV;
  obs = model.obs;
  K = model.K;
  alpha = model.alpha;
  C = model.C;
  Q = model.Q;
  current_size  = model.current_size ;


% first compute simple upate quantities
    k_t1 = kernel(x,BV,sigma)';   % pg 9, eqs 30-31   
    noise_x = noise + k_t1'*C*k_t1 + 1;
    q_t1 = (y - k_t1'*alpha)/(noise_x + noise);
    r_t1 = -1/(noise_x + noise);
    
    % compute residual projection update quantities 
    e_t1 = Q*k_t1; %residual vector pg 6, eq 16
    gamma_t1 = double(1-k_t1'*e_t1); %novelty of new point w.r.t RKHS: pg 7, eq 23
    eta_t1 = 1/(1+gamma_t1*r_t1); %numerical stability parameter    
    
    if gamma_t1 < tol
      % in this case, addition of point to basis doesn't help much, so 
      % don't add it, and compute update quantities in terms of old vectors
      % note that data, obs and gram matrix inverse not updated      
      s_t1 = C*k_t1 + e_t1;                  %pg 5, eqs 9, but modified
      alpha = alpha + q_t1*eta_t1*s_t1;
      C = C + r_t1*eta_t1*(s_t1*s_t1');                  
    else
      % in this case, you need to add the points
      current_size = current_size + 1;      
      
      %in this case, you can simply add the points    
      s_t1 = [C*k_t1; 1];    
      alpha = [alpha; 0] + q_t1*s_t1;
      C = [C zeros(current_size-1,1); zeros(1,current_size)] + r_t1*(s_t1*s_t1'); 
    
      % update basis vectors and observations
      BV = [BV x];
      obs = [obs; y];
    
      % update Gram matrix and inverse      
      K = [K k_t1; k_t1' 1]; 
      Q = inv(K);        
      
      if current_size <= max_points
        %do nothing   
      else
        % now you must delete one of the basis vectors; follow figure 3.3
        % first, compute which vector is least informative (1), pg 8, eq 27
        scores = zeros(1,current_size);        
        for i=1:current_size
         scores(i) = abs(alpha(i))/Q(i,i);   
        end
        
        %find index of minimum vector
        [val index] = min(scores);
        
        %now begin update given in (1), pg 8, eq 25
        
        %first compute scalar parameters
        a_s = alpha(index);
        c_s = C(index,index);
        q_s = Q(index,index);
        
        %compute vector parameters
        C_s = C(:,index);
        C_s(index) = [];
        Q_s = Q(:,index);
        Q_s(index) = [];

        %shrink matrices
        alpha(index) = [];
        C(:,index)   = [];
        C(index,:)   = [];
        Q(:,index)   = [];
        Q(index,:)   = [];
        K(:,index)   = [];
        K(index,:)   = [];
        
        %finally, compute updates
        alpha = alpha - (a_s/q_s)*(Q_s);
        C = C + (c_s/(q_s^2))*(Q_s*Q_s') - (1/q_s)*(Q_s*C_s' + C_s*Q_s');
        Q = Q - (1/q_s)*(Q_s*Q_s');
        
        current_size = current_size - 1;
        BV(:,index) = [];
        obs(index) = [];
      end    
      
%       c = (1/gamma_t1);
%       mc = -c;
%       Q = [(Q + c*(e_t1*e_t1')) mc*e_t1; mc*e_t1' c];          e 
     
    end 

    model_1 = model;

    model_1.BV = BV;
    model_1.obs = obs;
    model_1.current_size= current_size;
    model_1.Q = Q;
    model_1.C = C;
    model_1.K = K;
    model_1.alpha = alpha;




end

function v =  kernel(x,y,sigma)


  if(length(sigma) == 1) %same sigma
            d=x'*y;
            dx = sum(x.^2,1);
            dy = sum(y.^2,1);
            val = repmat(dx',1,length(dy)) + repmat(dy,length(dx),1) - 2*d;
            v = exp(-val./(2*sigma^2));
        else
            isigma = inv(diag(sigma.^2));
            d =  (x'*isigma)*y;
            dx = sum((x'*isigma)'.*x,1);
            dy = sum((y'*isigma)'.*y,1);
            val = repmat(dx',1,length(dy)) + repmat(dy,length(dx),1) - 2*d;
            v = exp(-val./2);
  end
  end