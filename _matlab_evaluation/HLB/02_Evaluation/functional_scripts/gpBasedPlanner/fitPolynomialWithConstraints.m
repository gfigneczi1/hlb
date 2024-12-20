function yFine = fitPolynomialWithConstraints(x,y, xFine, m)

    x0 = x(1);
    y0 = y(1);
    % 'C' is the Vandermonde matrix for 'x'
    V(:,m+1) = ones(size(x,1),1,class(x));
    for j = m:-1:1
       V(:,j) = x(:,1).*V(:,j+1);
    end
    C = V;
    % 'd' is the vector of target values, 'y'.
    d = y(:,1);
    %%
    % There are no inequality constraints in this case, i.e.,
    A = [];
    b = [];
    %%
    % We use linear equality constraints to force the curve to hit the required point. In
    % this case, 'Aeq' is the Vandermoonde matrix for 'x0'
    Aeq = x0.^(m:-1:0);
    % and 'beq' is the value the curve should take at that point
    beq = y0;
    %%
    p = lsqlin( C, d, A, b, Aeq, beq );
    %%
    % We can then use POLYVAL to evaluate the fitted curve
    yFine = polyval( p, xFine );

end
