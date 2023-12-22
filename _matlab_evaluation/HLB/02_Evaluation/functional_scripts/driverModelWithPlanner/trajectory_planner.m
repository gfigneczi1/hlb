function [traj, segments] = trajectory_planner(anchorPoints, indeces, ref, ver)
    %TRAJECTORY_PLANNER Summary of this function goes here
    %   N: number of segments (possible 1-4)
    %   c0Vector: Nx1 vector of c0 (start curvature) values
    %   c1Vector: Nx1 vector of c1 (curvature gradient) values
    %   Lvector: Nx1 vector of L-s (segment lengths)
    %   linkPoints: Nx2 matrix of linkPoints in the planning frame
    %   planning frame: ego frame at the time of planning
    %       Note: planning frame can be reinitalized in different ways (e.g. in
    %       all cycles iteratively or occasionally due to some external
    %       triggers
    % INPUT:
    %   corridor: is a clothoid of the mid-lane received from the GPS data in
    %   the form of 2d array of points (X-Y in the planning frame) m x 1, where
    %   m is the number of points
    %   ver is a development variable which means: ver = 0: clothoid, ver = 1
    %   3pcs of 5th order polynomials, ver =2: 3rd order polynomial for all
    %   node points ver = 3: cubic spline
    N = length(anchorPoints) / 3; % number of segments
    traj = []; segments = [];

    for i=1:N-1
        segmentStart = anchorPoints(i:N:end);
        segmentStop = anchorPoints(i+1:N:end);
        segment = segment_planner(segmentStart, segmentStop, ver);
        if (isempty(segment))
            segments = segment;
        else
            segments = [segments; segment];
        end
        switch ver
            case 0
                if (i==1)
                    % in case of the first segment, the clothoid fitting
                    % contains different number of points
                    XY = pointsOnClothoid(segment(1),segment(2), segment(3), ...
                     segment(4), segment(5), segment(6), indeces(i+1));
                    traj_temp = [ref(indeces(i):indeces(i+1)-1,1) spline(XY(1,1:indeces(i+1)),XY(2,1:indeces(i+1)),ref(1:indeces(i+1)-1,1))];
                elseif (i==N-1)
                    % last segment
                    XY = pointsOnClothoid(segment(1),segment(2), segment(3), ...
                     segment(4), segment(5), segment(6), indeces(i+1)-indeces(i));
                    traj_temp = [ref(indeces(i):indeces(i+1),1) spline(XY(1,1:indeces(i+1)-indeces(i)),XY(2,1:indeces(i+1)-indeces(i)),ref(indeces(i):indeces(i+1),1))];
                else
                    XY = pointsOnClothoid(segment(1),segment(2), segment(3), ...
                     segment(4), segment(5), segment(6), indeces(i+1)-indeces(i));
                    traj_temp = [ref(indeces(i):indeces(i+1)-1,1) spline(XY(1,1:indeces(i+1)-indeces(i)),XY(2,1:indeces(i+1)-indeces(i)),ref(indeces(i):indeces(i+1)-1,1))];
                end
                
                clear XY;
                if (isempty(traj))
                    traj = traj_temp;
                else
                    traj = [traj; traj_temp];
                end
        end
    end

%     initPose = anchorPoints(1:4:end);
%     mid1Pose = anchorPoints(2:4:end);
%     mid2Pose = anchorPoints(3:4:end);
%     endPose = anchorPoints(4:4:end);
%     
%     segment1 = segment_planner(initPose, mid1Pose, ver);
%     segment2 = segment_planner(mid1Pose, mid2Pose, ver);
%     segment3 = segment_planner(mid2Pose, endPose, ver);
%     segments = [segment1; segment2; segment3];
% 
%     switch ver
%         case 0
%             XY = pointsOnClothoid(segment1(1),segment1(2), segment1(3), ...
%                  segment1(4), segment1(5), segment1(6), indeces(2));
%             % p = polyfit(XY(1,:),XY(2,:),5);
%             % p(6) = XY(2,1);
%             traj_temp = [ref(1:indeces(2)-1,1) spline(XY(1,1:indeces(2)),XY(2,1:indeces(2)),ref(1:indeces(2)-1,1))]; %p(6)+p(5)*ref(indeces(2):indeces(3)-1,1)+p(4)*ref(indeces(2):indeces(3)-1,1).^2+p(3)*ref(indeces(2):indeces(3)-1,1).^3+p(2)*ref(indeces(2):indeces(3)-1,1).^4+p(1)*ref(indeces(2):indeces(3)-1,1).^5];
% 
%             clear XY;
%             traj = traj_temp;
%             XY = pointsOnClothoid(segment2(1),segment2(2), segment2(3), ...
%                  segment2(4), segment2(5), segment2(6), indeces(3)-indeces(2));
%              %p = polyfit(XY(1,:),XY(2,:),5);
%              traj_temp = [ref(indeces(2):indeces(3)-1,1) spline(XY(1,1:indeces(3)-indeces(2)),XY(2,1:indeces(3)-indeces(2)),ref(indeces(2):indeces(3)-1,1) )]; %p(6)+p(5)*ref(indeces(2):indeces(3)-1,1)+p(4)*ref(indeces(2):indeces(3)-1,1).^2+p(3)*ref(indeces(2):indeces(3)-1,1).^3+p(2)*ref(indeces(2):indeces(3)-1,1).^4+p(1)*ref(indeces(2):indeces(3)-1,1).^5];
%             %traj_temp(:,2) = traj_temp(:,2)  - (traj_temp(1,2)-XY(2,1)); 
%              traj = [traj; traj_temp];
%             clear XY;
% 
%             XY = pointsOnClothoid(segment3(1),segment3(2), segment3(3), ...
%                  segment3(4), segment3(5), segment3(6), indeces(4)-indeces(3));
%              %p = polyfit(XY(1,:),XY(2,:),5);
% 
%              traj_temp =[ref(indeces(3):indeces(4),1) spline(XY(1,1:indeces(4)-indeces(3)), XY(2,1:indeces(4)-indeces(3)), ref(indeces(3):indeces(4),1))]; % p(6)+p(5)*ref(indeces(3):indeces(4),1)+p(4)*ref(indeces(3):indeces(4),1).^2+p(3)*ref(indeces(3):indeces(4),1).^3+p(2)*ref(indeces(3):indeces(4),1).^4+p(1)*ref(indeces(3):indeces(4),1).^5];
%              %traj_temp(:,2) = traj_temp(:,2)  - (traj_temp(1,2)-XY(2,1)); 
%              traj = [traj; traj_temp];
% 
%              clear XY;
%         case 1
% %             theta = diff(ref(:,2))./diff(ref(:,1));
% %             kappa = diff([theta;theta(end)])./diff(ref(:,1));
% %             kappa = [kappa; kappa(end)];
% %             kappa = movmean(kappa,50);
% %             kappa0 = kappa(1);
% %             kappa1 = kappa(indeces(2));
% %             kappa2 = kappa(indeces(3));
% %             kappa3 = kappa(indeces(4));
% %             x = fit5thOrderPolynomial(initPose,mid1Pose,kappa0,kappa1);
% %             traj_temp = pointsOnPolynomial(x,ref(1:indeces(2)-1,1),initPose(1:2),mid1Pose(1:2));
% %             traj = traj_temp;
% %             x = fit5thOrderPolynomial(mid1Pose,mid2Pose,kappa1,kappa2);
% %             traj_temp = pointsOnPolynomial(x,ref(indeces(2):indeces(3)-1,1),mid1Pose(1:2),mid2Pose(1:2));
% %             traj = [traj; traj_temp];
% %             x = fit5thOrderPolynomial(mid2Pose,endPose,kappa2,kappa3);
% %             traj_temp = pointsOnPolynomial(x,ref(indeces(3):indeces(4),1),mid2Pose(1:2),endPose(1:2));
% %             traj = [traj; traj_temp];
%             x = fit5thOrderPolynomial(initPose,mid1Pose,mid2Pose,endPose);
%             traj = pointsOnPolynomial(x,ref(:,1),initPose(1:2),endPose(1:2));
%         case 2
%             x = polyfit(anchorPoints(1:4),anchorPoints(5:8),3);
%             x = flip(x);
%             x = [x'; 0; 0];
%             traj = pointsOnPolynomial(x,ref(:,1),initPose(1:2),endPose(1:2));
%         case 3
%             traj = [ref(1:indeces(4),1) spline(anchorPoints(1:4),anchorPoints(5:8),ref(1:indeces(4),1))];
% 
%     end
end

% function [x] = fit5thOrderPolynomial(initPose,endPose,kappa0,kappa1)
%     x0 = initPose(1); x1 = endPose(1);
%     y0 = initPose(2); y1 = endPose(2);
%     m0 = tan(initPose(3)); m1 = tan(endPose(3));
%     A = [1 x0 x0^2 x0^3 x0^4 x0^5; ...
%         0 1 2*x0 3*x0^2 4*x0^3 5*x0^4; ...
%         0 0 2 6*x0 12*x0^2 20*x0^3; ...
%         1 x1 x1^2 x1^3 x1^4 x1^5; ...
%         0 1 2*x1 3*x1^2 4*x1^3 5*x1^4; ...
%         0 0 2 6*x1 12*x1^2 20*x1^3];
%     b = [y0; m0; kappa0; y1; m1; kappa1];
%     %x = inv(A)*b;
%     x = A\b;
% end

function [x] = fit5thOrderPolynomial(initPose, mid1Pose, mid2Pose, endPose)
    x0 = initPose(1); x1 = mid1Pose(1); x2 = mid2Pose(1); x3 = endPose(1);
    y0 = initPose(2); y1 = mid1Pose(2); y2 = mid2Pose(2); y3 = endPose(2);
    m0 = tan(initPose(3)); m1 = tan(endPose(3));
    A = [1 x0 x0^2 x0^3 x0^4 x0^5; ...
        1 x1 x1^2 x1^3 x1^4 x1^5; ...
        1 x2 x2^2 x2^3 x2^4 x2^5; ...
        1 x3 x3^2 x3^3 x3^4 x3^5; ...
        0 1 2*x0 3*x0^2 4*x0^3 5*x0^4; ...
        0 1 2*x3 3*x3^2 4*x3^3 5*x3^4];
    b = [y0; y1; y2; y3; m0; m1];
    %x = inv(A)*b;
    x = A\b;
end

function XY = pointsOnPolynomial (coefficients, ref, initPoint, endPoint)
xpoints = ref(:,1);
initIdx = find(xpoints > initPoint(1),1);
xpoints = [initPoint(1);xpoints(initIdx:end)];
stopIdx = find (xpoints > endPoint(1),1);
if (~isempty(stopIdx))
    xpoints = [xpoints(1:stopIdx-1); endPoint(1)];
else
    xpoints = [xpoints; endPoint(1)];
end
ypoints = coefficients(1) + coefficients(2)*xpoints + coefficients(3)*xpoints.^2 + coefficients(4)*xpoints.^3 + coefficients(5)*xpoints.^4 + coefficients(6)*xpoints.^5;
XY = [xpoints(1:end-1) ypoints(1:end-1)];


end

function varargout = pointsOnClothoid( varargin )

  if nargin == 2
    if isstruct(varargin{1})
      x0     = varargin{1}.x0 ;
      y0     = varargin{1}.y0 ;
      theta0 = varargin{1}.theta0 ;
      kappa  = varargin{1}.kappa ;
      dkappa = varargin{1}.dkappa ;
      L      = varargin{1}.L ;
      npts   = varargin{2} ;
      tvec   = [0:L/(npts-1):L] ;
    else
      error('expexted struct as first arument') ;
    end
  elseif nargin == 6 || nargin == 7
    x0     = varargin{1} ;
    y0     = varargin{2} ;
    theta0 = varargin{3} ;
    kappa  = varargin{4} ;
    dkappa = varargin{5} ;
    L      = varargin{6} ;
    if nargin == 7
      npts = varargin{7} ;
      tvec = linspace(0,L,npts); %[0:L/(npts-1):L] ;
    else
      tvec = L ;
    end
  else
    error('expexted 2,6 or 7 input arguments') ;
  end

  X = zeros(1,25) ;
  Y = zeros(1,25) ;
  for i=1:length(tvec)
      t = tvec(i);
    [C,S] = GeneralizedFresnelCS( 1, dkappa*t^2, kappa*t, theta0 ) ;
    X(1,i) = x0 + t*C;
    Y(1,i) = y0 + t*S;
  end
  if nargout > 1
    varargout{1} = X ;
    varargout{2} = Y ;
  else
    varargout{1} = [X ; Y] ;
  end
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

