function phi = get_phi(gpmodel,newstate,flag)

% compute Kxm
[n, D] = size(gpmodel.inputs);    % number of examples and dimension of inputs
E = size(gpmodel.targets,2);         % number of examples and number of outputs
input = gpmodel.inputs; targets = gpmodel.targets;
%-----hyperparameter modified-------%
% X = gpmodel.hyp; 
X(1:D,:) = gpmodel.hyp(1:D,:)*(-0.5);
X(D+1:D+2,:) = gpmodel.hyp(D+1:end,:)*(0.5);

%-----hyperparameter modified-------%
[np pD pE] = size(gpmodel.induce);     % number of pseudo inputs per dimension
pinput = gpmodel.induce;     

if flag == 1 % 计算第一维

    pinp = bsxfun(@rdivide,pinput(:,:,1),exp(X(1:D,1)'));
    pnew = bsxfun(@rdivide,newstate,exp(X(1:D,1)'));
    phi =  exp(2*X(D+1,1)-maha(pnew,pinp)/2); 
elseif flag == 2
    pinp = bsxfun(@rdivide,pinput(:,:,2),exp(X(1:D,2)'));
    pnew = bsxfun(@rdivide,newstate,exp(X(1:D,2)'));
    phi =  exp(2*X(D+1,2)-maha(pnew,pinp)/2);
elseif flag == 3
    pinp = bsxfun(@rdivide,pinput(:,:,3),exp(X(1:D,3)'));
    pnew = bsxfun(@rdivide,newstate,exp(X(1:D,3)'));
    phi =  exp(2*X(D+1,3)-maha(pnew,pinp)/2); 
end