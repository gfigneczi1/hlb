function v =  kernel_o(x,y,sigma)


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