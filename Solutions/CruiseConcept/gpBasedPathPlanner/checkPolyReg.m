X = [rand(1,100); rand(1,100)]';
c = [0.1 0.15 0.2 0.25 0.3 0.35]';
y = c(1:2)'*X'+c(3:4)'*(X.^2)'+c(5:6)'*(X.^3)';

c_ = polynomialRegression(X,y',3);
y_ = c_(1:2)'*X'+c_(3:4)'*(X.^2)'+c_(5:6)'*(X.^3)';

function coefficients = polynomialRegression(X,y, p)
% formula from here: https://math.stackexchange.com/questions/3155866/multivariate-quadratic-regression
% mathematically: (A'A)^(-1)*A'*y
% where: A[i,:] = [1 xi1 xi2 ... xin xi1^2 xi2^2 ... xin^2 ...xin^p] where
% p is the order of the polynomial
for i = 1:size(X,1)
    % loop through input samples
    x = X(i,:);
    for j=1:p
        A(i,(j-1)*size(X,2)+1:j*size(X,2)) = x.^(j);
    end
end
coefficients = inv(A'*A)*A'*y;
end