%% fitc.m
% *Summary:* Compute the FITC negative log marginal likelihood and its
% derivatives with  respect to the inducing inputs (we don't compute the
% derivatives with respect to the GP hyper-parameters)
%
%   function [nml dnml] = fitc(induce, gpmodel)
%
% *Input arguments:*
%
%   induce          matrix of inducing inputs                       [M x D x uE]
%                   M: number of inducing inputs
%                   E: either 1 (inducing inputs are shared across target dim.)
%                      or     E (different inducing inputs for each target dim.)
%   gpmodel         GP structure
%     .hyp          log-hyper-parameters                               [D+2 x E]
%     .inputs       training inputs                                    [N   x D]
%     .targets      training targets                                   [N   x E]
%     .noise (opt)  noise
%
%
% *Output arguments:*
%
%   nlml             negative log-marginal likelihood
%   dnlml            derivative of negative log-marginal likelihood wrt
%                    inducing inputs
%
% Adapted from Ed Snelson's SPGP code.
%
% Copyright (C) 2008-2013 by
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modified: 2013-05-21
%
%% High-Level Steps
%
% # Compute FITC marginal likelihood 
% # Compute corresponding gradients wrt the pseudo inputs
% # and the hyps

function [nlml dnlml] = fitc_mine(w, gpmodel)
%% Code

del = 1e-06;                       % jitter to make matrix better conditioned

[N,dim] = size(gpmodel.inputs); E = size(gpmodel.targets,2);
[M uD uE] = size(gpmodel.induce);
n = M;
for i = 1:uE
induce(:,:,i) = reshape(w(1:end-uD-2,i),M,uD);
end


if uD ~= dim || (uE~=1 && uE ~= E); error('Wrong size of inducing inputs'); end

nlml = 0; dfxb = zeros(M,dim); dnlml = zeros(M*dim+dim+2, E); % zero and allocate outputs

for j = 1:E
  if uE > 1; xb = induce(:,:,j); else xb = induce; end
  b = exp(w(end-dim-1:end-2,j)); c = exp(w(end-1,j)); sig = exp(w(end,j));

xb = xb.*repmat(sqrt(b)',n,1);
x = gpmodel.inputs.*repmat(sqrt(b)',N,1);
y = gpmodel.targets(:,j);                                  % training targets

Q = xb*xb';
Q = repmat(diag(Q),1,n) + repmat(diag(Q)',n,1) - 2*Q;
Q = c*exp(-0.5*Q) + del*eye(n);

K = -2*xb*x' + repmat(sum(x.*x,2)',n,1) + repmat(sum(xb.*xb,2),1,N);
K = c*exp(-0.5*K);

L = chol(Q)';
V = L\K;
ep = 1 + (c-sum(V.^2)')/sig;
K = K./repmat(sqrt(ep)',n,1);
V = V./repmat(sqrt(ep)',n,1); y = y./sqrt(ep);
Lm = chol(sig*eye(n) + V*V')';
invLmV = Lm\V;
bet = invLmV*y;

% Likelihood
nlml = nlml + sum(log(diag(Lm))) + (N-n)/2*log(sig) + ... 
      (y'*y - bet'*bet)/2/sig + sum(log(ep))/2 + 0.5*N*log(2*pi);

  
   if nargout == 2               % ... and if requested, its partial derivatives

  % precomputations
    Lt = L*Lm;
    B1 = Lt'\(invLmV);
    b1 = Lt'\bet;
    invLV = L'\V;
    invL = inv(L); invQ = invL'*invL; clear invL
    invLt = inv(Lt); invA = invLt'*invLt; clear invLt
    mu = ((Lm'\bet)'*V)';
    sumVsq = sum(V.^2)'; clear V
    bigsum = y.*(bet'*invLmV)'/sig - sum(invLmV.*invLmV)'/2 - (y.^2+mu.^2)/2/sig ...
             + 0.5;
    TT = invLV*(invLV'.*repmat(bigsum,1,n));

    % pseudo inputs and lengthscales
    for i = 1:dim
    % dnnQ = (repmat(xb(:,i),1,n)-repmat(xb(:,i)',n,1)).*Q;
    % dNnK = (repmat(x(:,i)',n,1)-repmat(xb(:,i),1,N)).*K;
    dnnQ = dist(xb(:,i),xb(:,i)).*Q;  % 20*20
    dNnK = dist(-xb(:,i),-x(:,i)).*K;

    epdot = -2/sig*dNnK.*invLV; epPmod = -sum(epdot)';

    dfxb(:,i) = - b1.*(dNnK*(y-mu)/sig + dnnQ*b1) ...
        + sum((invQ - invA*sig).*dnnQ,2) ...
        + epdot*bigsum - 2/sig*sum(dnnQ.*TT,2); 

    dfb(i,1) = (((y-mu)'.*(b1'*dNnK))/sig ...
               + (epPmod.*bigsum)')*x(:,i);

    dNnK = dNnK.*B1; % overwrite dNnK
    dfxb(:,i) = dfxb(:,i) + sum(dNnK,2);
    dfb(i,1) = dfb(i,1) - sum(dNnK,1)*x(:,i);

    dfxb(:,i) = dfxb(:,i)*sqrt(b(i));

    dfb(i,1) = dfb(i,1)/sqrt(b(i));
    dfb(i,1) = dfb(i,1) + dfxb(:,i)'*xb(:,i)/b(i);
    dfb(i,1) = dfb(i,1)*sqrt(b(i))/2;
    end

    % size
    epc = (c./ep - sumVsq - del*sum((invLV).^2)')/sig;

    dfc = (n + del*trace(invQ-sig*invA) ... 
         - sig*sum(sum(invA.*Q')))/2 ...
        - mu'*(y-mu)/sig + b1'*(Q-del*eye(n))*b1/2 ... 
          + epc'*bigsum;

    % noise
    dfsig = sum(bigsum./ep);

    dfw = [reshape(dfxb,n*dim,1);dfb;dfc;dfsig];
    %     
    dnlml(:,j) = dfw;
    
    end
end
% if 1 == uE; dnlml = sum(dnlml,3); end % combine derivatives if sharing inducing
