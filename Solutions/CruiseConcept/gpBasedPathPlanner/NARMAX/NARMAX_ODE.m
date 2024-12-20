%------------------------------------------------------------------
function [dx, y] = NARMAX_ODE(~,x,u,fcn_par,C,varargin) 
%NARMAX_ODE ODE function for refining the coefficients of a NARMAX model 
% using a grey-box approach.
%
% Ther NARMAX model uses an additive moving average term.

% Copyright 2024 The MathWorks, Inc.

% If there are 2 moving average terms, the states are: 
% x1 := y(t-1)
% x2 := y(t-2)
% x3 := y_measured(t-1)
% x4 := y_measured(t-2)
% x5 := u(t-1)
% x6 := u(t-2)
% x7 := u(t-3)

% fcn_par -> coefficients used in f(.)
% C -> free coefficients of the moving average polynomial C(q)
% sys -> Nonlinear ARX model representing the function f(.)

% Separate out the parameters of the nonlinear sigmoidal function 
% and those corresponding to the moving average term
sys = varargin{1}{1};
nc = varargin{1}{2};
% sys = setpvec(sys, fcn_par); % not needed if fcn_par is fixed in the example
NL = sys.OutputFcn;
regressor_vector = x([nc+(1:2), end-2:end])';
y_NL = evaluate(NL,regressor_vector);
y_MA = C'*(x(nc+(1:nc))-x(1:nc));
y = y_NL + y_MA;
 
dx1 = y;
for k=1:nc-1
   dx1 = [dx1; x(k)]; %#ok<*AGROW>
end
dx1(end+1) = u(1);
for k=nc+(1:nc-1)
   dx1 = [dx1; x(k)];
end
dx2 = [u(2); x(end-2); x(end-1)];
dx = [dx1;dx2];

end
%------------------------------------------------------------------