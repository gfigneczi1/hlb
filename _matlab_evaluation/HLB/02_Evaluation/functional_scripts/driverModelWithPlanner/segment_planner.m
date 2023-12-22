function [segmentParams, coefficients] = segment_planner(initPose, endPose, ver)
%SEGMENT_PLANNER Summary of this function goes here
%   Detailed explanation goes here
switch ver
    case 0
        [k, dk, L] = buildClothoid (initPose(1), initPose(2), initPose(3), endPose(1), endPose(2), endPose(3));
        segmentParams = [initPose(1) initPose(2) initPose(3) k dk L];
    otherwise
        coefficients = buildPolynomial(initPose(1), initPose(2), endPose(1), endPose(2),initPose(3), endPose(3));
        segmentParams = [initPose(1) initPose(2) initPose(3) coefficients(3)*2 coefficients(4)*6 0];
end
end

function [coefficients] = buildPolynomial (x0, y0, x1, y1, theta0, theta1)

A = [1 x0 x0^2 x0^3; 0 1 2*x0 3*x0^2; 1 x1 x1^2 x1^3; 0 1 2*x1 3*x1^2];
b = [y0; tan(theta0); y1; tan(theta1)];
%x = inv(A)*b;
x = A\b;
coefficients(1) = x(1);
coefficients(2) = x(2);
coefficients(3) = x(3);
coefficients(4) = x(4);

end

function [ k, dk, L, iter, varargout ] = buildClothoid( x0, y0, theta0, x1, y1, theta1 )

  dx  = x1 - x0 ;
  dy  = y1 - y0 ;
  r   = sqrt( dx^2 + dy^2 ) ;
  phi = atan2( dy, dx ) ;

  phi0  = normalizeAngle(theta0 - phi) ;
  phi1  = normalizeAngle(theta1 - phi) ;
  delta = phi1 - phi0 ;

  % initial point
  Aguess = guessA( phi0, phi1 ) ;

  % Newton iteration
  [A,iter] = findA( Aguess, delta, phi0, 1e-12 ) ;

  % final operation
  [h,g] = GeneralizedFresnelCS( 1, 2*A, delta-A, phi0 ) ;
  L = r/h ;

  if L > 0
    k  = (delta - A)/L ;
    dk = 2*A/L^2 ;
  else
    error('negative length') ;
  end

  if nargout == 10

    [X,Y] = GeneralizedFresnelCS( 3, 2*A, delta-A, theta0 ) ;
    
    if true
      alpha = X(1)*X(2) + Y(1)*Y(2) ;
      beta  = X(1)*X(3) + Y(1)*Y(3) ;
      gamma = X(1)^2+Y(1)^2 ;
      tx    = X(2)-X(3) ;
      ty    = Y(2)-Y(3) ;
      txy   = L*(X(2)*Y(3)-X(3)*Y(2)) ;
      omega = L*(Y(1)*tx-X(1)*ty) - txy ;
      delta = X(1)*tx + Y(1)*ty ;

      L_1  = omega/delta ; % L_0
      L_2  = txy/delta ; % L_1

      delta = delta * L ;
      k_1  = (beta-gamma-k*omega)/delta ; % k_0
      k_2  = -(beta+k*txy)/delta ; % k_1

      delta = delta * L/2 ;
      dk_1 = (gamma-alpha-dk*omega*L)/delta ; % dk_0    
      dk_2 = (alpha-dk*txy*L)/delta ; % dk_1
    else
      dkL = dk*L ;
      M = [ X(1) - L*(dkL*Y(3)+k*Y(2)), -L*Y(2), -L*Y(3) ; ...
            Y(1) + L*(dkL*X(3)+k*X(2)),  L*X(2),  L*X(3) ; ...
            dkL+k,                            1,       1 ] ;
      tmp = M\[L*Y(1);-L*X(1);-1] ;
      L_1  = tmp(1) ;
      k_1  = tmp(2)/L ;
      dk_1 = 2*tmp(3)/L^2 ;
      tmp = M\[0;0;1] ;
      L_2  = tmp(1) ;
      k_2  = tmp(2)/L ;
      dk_2 = 2*tmp(3)/L^2 ;
    end
    
    varargout{1} = k_1 ; % k_0
    varargout{2} = dk_1 ; % dk_0
    varargout{3} = L_1 ; % L_0
    
    varargout{4} = k_2 ; % k_1
    varargout{5} = dk_2 ; % dk_1
    varargout{6} = L_2 ; % L_1

  elseif nargout > 4
    error('expected <= 4 or 10 output argument') ;
  end
  
end

%=========================================================================%
%  normalizeAngle:  normalize angle in the range [-pi,pi]                 %
%=========================================================================%
function phi = normalizeAngle( phi_in )
  phi = phi_in ;
  while ( phi > pi )
    phi = phi - 2*pi ;
  end
  while ( phi < -pi )
    phi = phi + 2*pi ;
  end
end

%=========================================================================%
%  findA:  Find a zero of function g(A) defined as                        %
%  g(A) = \int_0^1 \sin( A*t^2+(delta-A)*t+phi0 ) dt                      %
%                                                                         %
%  USAGE:  A = findA( Aguess, delta, phi0, tol );                         %
%                                                                         %
%  Given an initial guess Aguess find the closest zero of equation g(A)   %
%                                                                         %
%  On input:                                                              %
%    Aguess      = initial guess.                                         %
%    delta, phi0 = Angles used in the clothoid fitting problem.           %
%    tol         = Tolerance for stopping criterium of Newton iteration.  %
%                                                                         %
%  On output:                                                             %
%    A           = the zero of function g(A) closest to Aguess.           %
%    iter        = iteration performed                                    %
%                                                                         %
%=========================================================================%
function [A,iter] = findA( Aguess, delta, phi0, tol )
  A = Aguess ;
  for iter=1:100
    [intC,intS] = GeneralizedFresnelCS( 3, 2*A, delta-A, phi0 ) ;
    f  = intS(1) ;
    df = intC(3)-intC(2) ;
    A  = A - f/df ;
    if abs(f) < tol
      break ;
    end
  end
  if abs(f) > tol*10
    fprintf( 1, 'Newton iteration fails, f = %g\n', f ) ;
    fprintf( 1, 'Aguess = %g, A = %g, delta = %g , phi0 = %g\n', Aguess, A, delta, phi0 ) ;
  end
end

%=========================================================================%
%  guessA:  Find guess for zeros of function g(A)                         %
%                                                                         %
%  USAGE:  A = guessA( phi0, phi1 );                                      %
%                                                                         %
%  On input:                                                              %
%       phi0, phi1 = Angles used in the clothoid fitting problem.         %
%                                                                         %
%  On output:                                                             %
%       A = an approximate zero of function g(A).                         %
%                                                                         %
%=========================================================================%
function A = guessA( phi0, phi1 )
  CF = [ 2.989696028701907, ...
         0.716228953608281, ...
        -0.458969738821509, ...
        -0.502821153340377, ...
         0.261062141752652, ...
        -0.045854475238709 ] ;
  X  = phi0/pi ;
  Y  = phi1/pi ;
  xy = X*Y ;
  A  = (phi0+phi1) * ( CF(1) + xy * ( CF(2) + xy * CF(3) ) + ...
                       (CF(4) + xy * CF(5)) * (X^2+Y^2) + ...
                        CF(6) * (X^4+Y^4) ) ;
end

%=============================================================================%
%  intXY:  Compute Fresnel sine and cosine integrals momenta                  %
%                                                                             %
%  USAGE: [X,Y] = GeneralizedFresnelCS( nk, a, b, c ) ;                       %
%                                                                             %
%  Integrals are defined as:                                                  %
%                                                                             %
%  X_k(a,b,c) = \int_0^1 t^k * cos( (a/2)*t^2 + b*t + c ) dt                  %
%  Y_k(a,b,c) = \int_0^1 t^k * sin( (a/2)*t^2 + b*t + c ) dt                  %
%                                                                             %
%  On input:                                                                  %
%                                                                             %
%       nk      = number of momentae to be computed                           %
%       a, b, c = the parameters of the integrals                             %
%                                                                             %
%  On output:                                                                 %
%                                                                             %
%       X = vector with Fresnel cosine momenta [X_0,X_1,...,X_{nk-1}]         %
%       Y = vector with Fresnel sine momenta   [Y_0,Y_1,...,Y_{nk-1}]         %
%                                                                             %
%=============================================================================%
%                                                                             %
%  Autors: Enrico Bertolazzi and Marco Frego                                  %
%          Department of Industrial Engineering                               %
%          University of Trento                                               %
%          enrico.bertolazzi@unitn.it                                         %
%          m.fregox@gmail.com                                                 %
%                                                                             %
%=============================================================================%
function [X,Y] = GeneralizedFresnelCS( nk, a, b, c )

  epsi = 1e-2 ; % best thresold

  if abs(a) < epsi % case `a` small
    [X,Y] = evalXYaSmall( nk, a, b, 3 ) ;
  else
    [X,Y] = evalXYaLarge( nk, a, b ) ;
  end

  cc = cos(c) ;
  ss = sin(c) ;

  for k=1:nk
    xx = X(k) ;
    yy = Y(k) ;
    X(k) = xx*cc-yy*ss ;
    Y(k) = xx*ss+yy*cc ;
  end

end
%
%
%
function [C,S] = FresnelCSk( nk, t )
  C = zeros(nk,1) ;
  S = zeros(nk,1) ;
  [C(1),S(1)] = FresnelCS(t) ;
  if nk > 1
    tt   = pi/2*t^2 ;
    ss   = sin(tt) ;
    cc   = cos(tt) ;
    C(2) = ss/pi ;
    S(2) = (1-cc)/pi ;
    if nk > 2
      C(3) = (t*ss-S(1))/pi ;
      S(3) = (C(1)-t*cc)/pi ;
    end
  end
end
%
%
%
function [X,Y] = evalXYaLarge( nk, a, b )
  X       = zeros(nk,1) ;
  Y       = zeros(nk,1) ;
  s       = sign(a) ;
  z       = sqrt(abs(a)/pi) ;
  ell     = s*b/sqrt(abs(a)*pi) ;
  g       = -0.5*s*b^2/abs(a) ;
  [Cl,Sl] = FresnelCSk(nk,ell) ;
  [Cz,Sz] = FresnelCSk(nk,ell+z) ;
  dC      = Cz - Cl ;
  dS      = Sz - Sl ;
  cg      = cos(g)/z ;
  sg      = sin(g)/z ;
  X(1)    = cg * dC(1) - s * sg * dS(1) ;
  Y(1)    = sg * dC(1) + s * cg * dS(1) ;
  if nk > 1
    cg   = cg/z ;
    sg   = sg/z ;
    DC   = dC(2)-ell*dC(1);
    DS   = dS(2)-ell*dS(1);
    X(2) = cg * DC - s * sg * DS ;
    Y(2) = sg * DC + s * cg * DS ;
    if nk > 2
      cg   = cg/z ;
      sg   = sg/z ;
      DC   = dC(3)+ell*(ell*dC(1)-2*dC(2)) ;
      DS   = dS(3)+ell*(ell*dS(1)-2*dS(2)) ;
      X(3) = cg * DC - s * sg * DS ;
      Y(3) = sg * DC + s * cg * DS ;
    end
  end
end
%
%
%
function res = rLommel_MATLAB(mu, nu, z)
  a   = mu-nu+1 ;
  b   = mu+nu+1 ;
  res = eval(hypergeom(1,[a/2+1,b/2+1],-z^2/4)/(a*b)) ;
end
%
%
%
function res = rLommel(mu,nu,b)
  tmp = 1/((mu+nu+1)*(mu-nu+1)) ;
  res = tmp ;
  for n=1:100
    tmp = tmp * (-b/(2*n+mu-nu+1)) * (b/(2*n+mu+nu+1)) ;
    res = res + tmp ;
    if abs(tmp) < abs(res) * 1e-50
      break ;
    end;
  end
end
%
%
%
function [X,Y] = evalXYazero( nk, b )
  X  = zeros(nk,1) ;
  Y  = zeros(nk,1) ;
  sb = sin(b);
  cb = cos(b);
  b2 = b*b ;
  if abs(b) < 1e-3
    X(1) = 1-(b2/6)*(1-(b2/20)*(1-(b2/42))) ;
    Y(1) = (b/2)*(1-(b2/12)*(1-(b2/30))) ;
  else
    X(1) = sb/b ;
    Y(1) = (1-cb)/b ;
  end
  % use recurrence in the stable part
  m = min( [ max([1,floor(2*b)]), nk] ) ;
  for k=1:m-1
    X(k+1) = (sb-k*Y(k))/b ;
    Y(k+1) = (k*X(k)-cb)/b ;
  end
  % use Lommel for the unstable part
  if m < nk
    A   = b*sb ;
    D   = sb-b*cb ;
    B   = b*D ;
    C   = -b2*sb ;
    rLa = rLommel(m+1/2,3/2,b) ;
    rLd = rLommel(m+1/2,1/2,b) ;
    for k=m:nk-1
      rLb    = rLommel(k+3/2,1/2,b) ;
      rLc    = rLommel(k+3/2,3/2,b) ;
      X(k+1) = ( k*A*rLa + B*rLb + cb )/(1+k) ;
      Y(k+1) = ( C*rLc + sb ) / (2+k) + D*rLd ;
      rLa = rLc ;
      rLd = rLb ;
    end
  end
end
%
%
%
function [X,Y] = evalXYaSmall( nk, a, b, p )
  [X0,Y0] = evalXYazero( nk + 4*p + 2, b ) ;

  X    = zeros(nk,1) ;
  Y    = zeros(nk,1) ;
  tmpX = zeros(p+1,1) ;
  tmpY = zeros(p+1,1) ;

  for j=1:nk
    tmpX(1) = X0(j)-(a/2)*Y0(j+2) ;
    tmpY(1) = Y0(j)+(a/2)*X0(j+2) ;
    t  = 1 ;
    aa = -(a/2)^2 ;
    for n=1:p
      ii = 4*n+j ;
      t  = t*(aa/(2*n*(2*n-1))) ;
      bf = a/(4*n+2) ;
      tmpX(1+n) = t*(X0(ii)-bf*Y0(ii+2)) ;
      tmpY(1+n) = t*(Y0(ii)+bf*X0(ii+2)) ;
    end
    X(j) = sum(tmpX) ;
    Y(j) = sum(tmpY) ;
  end
end

function [FresnelC,FresnelS] = FresnelCS(y)
%******************************************************************************
% FresnelCS:  Compute Fresnel sine and cosine integrals                       %
%                                                                             %
% USAGE: [FresnelC,FresnelS] = FresnelCS(y) ;                                 %
%                                                                             %
%  Fresnel integral are defined as:                                           %
%                                                                             %
%  C(y) = \int_0^y cos( (pi/2)*t^2 ) dt                                       %
%  S(y) = \int_0^y sin( (pi/2)*t^2 ) dt                                       %
%                                                                             %
% The algorithm is described in:                                              %
%   Atlas for computing mathematical functions : an illustrated guide for     %
%   practitioners, with programs in C and Mathematica / William J. Thompson.  %
%   New York : Wiley, c1997.                                                  %
%                                                                             %
% The cose is a sligly modification of original code developed by             %
% Venkata Sivakanth Telasula (sivakanth.telasula@gmail.com)                   %
%                                                                             %
% On input:                                                                   %
%                                                                             %
%       y = argument of the function C(y) and S(y), it may be a vector        %
%                                                                             %
% On output:                                                                  %
%                                                                             %
%       FresnelC = The value(s) of C(y)                                       %
%       FresnelS = The value(s) of S(y)                                       %
%                                                                             %
%******************************************************************************

  fn = [ 0.49999988085884732562,   ...
         1.3511177791210715095,    ...
         1.3175407836168659241,    ...
         1.1861149300293854992,    ...
         0.7709627298888346769,    ...
         0.4173874338787963957,    ...
         0.19044202705272903923,   ...
         0.06655998896627697537,   ...
         0.022789258616785717418,  ...
         0.0040116689358507943804, ...
         0.0012192036851249883877 ] ;
  
  fd = [ 1.0,                      ...
         2.7022305772400260215,    ...
         4.2059268151438492767,    ...
         4.5221882840107715516,    ...
         3.7240352281630359588,    ...
         2.4589286254678152943,    ...
         1.3125491629443702962,    ...
         0.5997685720120932908,    ...
         0.20907680750378849485,   ...
         0.07159621634657901433,   ...
         0.012602969513793714191,  ...
         0.0038302423512931250065 ] ;
  
  gn = [ 0.50000014392706344801,    ...
         0.032346434925349128728,   ...
         0.17619325157863254363,    ...
         0.038606273170706486252,   ...
         0.023693692309257725361,   ...
         0.007092018516845033662,   ...
         0.0012492123212412087428,  ...
         0.00044023040894778468486, ...
        -8.80266827476172521e-6,    ...
        -1.4033554916580018648e-8,  ...
        2.3509221782155474353e-10 ] ;
  
  gd  = [ 1.0,                      ...
          2.0646987497019598937,    ...
          2.9109311766948031235,    ...
          2.6561936751333032911,    ...
          2.0195563983177268073,    ...
          1.1167891129189363902,    ...
          0.57267874755973172715,   ...
          0.19408481169593070798,   ...
          0.07634808341431248904,   ...
          0.011573247407207865977,  ...
          0.0044099273693067311209, ...
         -0.00009070958410429993314 ] ;

  FresnelC = zeros(size(y));
  FresnelS = zeros(size(y));
  
  for j = 1:length(y)
    x = abs(y(j)) ;
    if x < 1.0
      t = -((pi/2)*x.^2).^2;
      % Cosine integral series
      twofn   = 0.0 ;
      fact    = 1.0 ;
      denterm = 1.0 ;
      numterm = 1.0 ;
      sum     = 1.0 ;
      ratio   = 10.0 ; 
 
      while ratio > eps
        twofn   = twofn + 2.0;
        fact    = fact*twofn*(twofn-1.0);
        denterm = denterm + 4.0;
        numterm = numterm*t;
        term    = numterm/(fact*denterm);
        sum     = sum+term;
        ratio   = abs(term/sum);
      end

      FresnelC(j) = x*sum;
 
      % Sine integral series
      twofn   = 1.0 ;
      fact    = 1.0 ;
      denterm = 3.0 ;
      numterm = 1.0 ;
      sum     = 1.0/3.0 ;
      ratio   = 10.0 ;
 
      while ratio > eps
        twofn   = twofn+2.0;
        fact    = fact*twofn*(twofn-1.0);
        denterm = denterm+4.0;
        numterm = numterm*t;
        term    = numterm/(fact*denterm);
        sum     = sum+term;
        ratio   = abs(term/sum);
      end
 
      FresnelS(j) = (pi/2)*sum.*x.^3;

    elseif x < 6.0

      % Rational approximation for f
      sumn = 0.0 ;
      sumd = fd(12) ;
      for k=11:-1:1
        sumn = fn(k)+x*sumn ;
        sumd = fd(k)+x*sumd ;
      end
      f = sumn/sumd;
      % Rational approximation for  g
      sumn = 0.0 ;
      sumd = gd(12) ;
      for k=11:-1:1
        sumn = gn(k)+x*sumn;
        sumd = gd(k)+x*sumd;
      end
      g    = sumn/sumd;
      U    = (pi/2)*x.^2 ;
      SinU = sin(U);
      CosU = cos(U);
      FresnelC(j) = 0.5+f*SinU-g*CosU;
      FresnelS(j) = 0.5-f*CosU-g*SinU;
    else
      % x >= 6; asymptotic expansions for  f  and  g
      t = -(pi*x.^2)^(-2.0) ;
      % Expansion for  f
      numterm = -1.0 ;
      term    =  1.0 ;
      sum     =  1.0 ;
      oldterm =  1.0 ;
      ratio   = 10.0 ;
      eps10   = 0.1*eps;

      while ratio > eps10
        numterm = numterm+4.0;
        term    = term*numterm*(numterm-2.0)*t;
        sum     = sum + term;
        absterm = abs(term);
        ratio   = abs(term/sum);
        if oldterm < absterm
          disp('\n\n !!In FresnelCS f not converged to eps');
          ratio = eps10;
        end
        oldterm = absterm;
      end

      f = sum/(pi*x);
      % Expansion for  g
      numterm = -1.0 ;
      term    =  1.0 ;
      sum     =  1.0 ;
      oldterm =  1.0 ; 
      ratio   = 10.0 ;
      eps10   = 0.1*eps ;

      while ratio > eps10
        numterm = numterm+ 4.0;
        term    = term*numterm*(numterm+2.0)*t;
        sum     = sum+term;
        absterm = abs(term);
        ratio   = abs(term/sum);
        if oldterm < absterm
          disp('\n\n!!In FresnelCS g not converged to eps');
          ratio = eps10;
        end
        oldterm = absterm;
      end
      g    = sum/((pi*x)^2*x) ;
      U    = (pi/2)*x.^2 ;
      SinU = sin(U);
      CosU = cos(U);
      FresnelC(j) = 0.5+f*SinU-g*CosU;
      FresnelS(j) = 0.5-f*CosU-g*SinU;
    end
    if y(j) < 0
      FresnelC(j) = -FresnelC(j) ;
      FresnelS(j) = -FresnelS(j) ;
    end
  end
end
% EOF
