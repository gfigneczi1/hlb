function [loss, gradients] = modelLoss(parameters, X, T)
    path = pathGeneration();
    
    loss = (path(:,1)-T(:,1)).^2+(path(:,2)-T(:,2)).^2;
    
    gradients = dlgradient(loss,parameters);
end

