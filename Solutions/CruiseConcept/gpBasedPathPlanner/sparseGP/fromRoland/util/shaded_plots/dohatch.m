function h = dohatch(x,y,angle,color,style,step,width)
%DOHATCH  Hatches a two-dimensional domain.
%   DOHATCH(x,y,ANG,COL,'Style',D,W) is similar to the FILL command
%   but fills the closed contour with hatch lines instead of uniform
%   color. Vectors x and y are the coordinates of the boundary line of
%   the domain to be hatched. Scalar ANG is the slope of the hatches
%   (in degrees). COL is the color of the hatching which could both be
%   a 1-by-3 vector ([red green blue] ), or a color specifier ('r','g',
%   'b', 'w','y','c','m'). 'Style' specifys the linestyle ('-','--','-.',
%   ':'). And also, D is the steps (distances between hatches), W the
%   linewidth (thickness) of the hatch lines (the last two in points).
%
%   H = DOHATCH(...) returns the handle of the hatching.
%
%   Note:
%   This is a simplified version of the function HATCH (File ID: #2075)
%   in MathWorks FileExchange website, running in manual style instead of
%   in automatic style, because the original version still has a few bugs
%   under some circumstances.
%
%   See also: FILL, LINE.

%   Revised by Kastin
%   June, 28, 2015

if nargin>7||nargin<2
    error('MATLAB:dohatch:WrongNumOfInputs',...
        'Wrong number of inputs.')
end

% Here we do not check the validation of inputs any more for clearity. 

% Defaults
angledflt = 45;        % Angle in degrees
colordflt = [1 1 1];   % Color
styledflt = '-';       % Linestyle
widthdflt = 1;         % Thickness of the lines
stepdflt = 10;         % Distance between hatches

if nargin<7, width = widthdflt; end
if nargin<6, step = stepdflt;   end
if nargin<5, style = styledflt; end
if nargin<4, color = colordflt; end
if nargin<3, angle = angledflt; end

angle = angle*pi/180;             % Degrees to radians
x=x(:).'; y=y(:).';
yi = find(~isnan(x)&~isnan(y));
x = x(yi); y = y(yi);             % Remove NaN's
x = [x x(1)]; y = [y y(1)];       % Close loop
ll = length(x);

% Transform the coordinates
oldu = get(gca,'units');
set(gca,'units','points')
sza = get(gca,'pos'); sza = sza(3:4);
xlim = get(gca,'xlim');
ylim = get(gca,'ylim');

islx = strcmp(get(gca,'xscale'),'log');
isly = strcmp(get(gca,'yscale'),'log');
if islx     % If log scale in x
	xlim = log10(xlim);
	x = log10(x);
end
if isly     % If log scale in y
	ylim = log10(ylim);
	y = log10(y);
end

xsc = sza(1)/(xlim(2)-xlim(1)+eps);
ysc = sza(2)/(ylim(2)-ylim(1)+eps);

ca = cos(angle); sa = sin(angle);
x0 = mean(x); y0 = mean(y);  % Central point
x = (x-x0)*xsc; y = (y-y0)*ysc;
yi = x*ca+y*sa;              % Rotation
y = -x*sa+y*ca;
x = yi;
y = y/step;    % Make steps equal to one

% Compute the coordinates of the hatch line
yi = ceil(y);
ll = length(y);
yd = [yi(2:ll)-yi(1:ll-1) 0];
dm = max(abs(yd));
fnd = find(yd);
lfnd = length(fnd);
A = sign(yd(fnd));
edm = ones(dm,1);
A = A(edm,:);
if size(A,1)>1, A = cumsum(A); end
fnd1 = find(abs(A)<=abs(yd(edm,fnd)));
A = A+yi(edm,fnd)-(A>0);
xy = (x(fnd+1)-x(fnd))./(y(fnd+1)-y(fnd));
xi = x(edm,fnd)+(A-y(edm,fnd)).*xy(edm,:);
yi = A(fnd1);
xi = xi(fnd1);

% Sorting points of the hatch line ........................
li = length(xi); 
xi0 = min(xi); xi1 = max(xi);
yi0 = min(yi); yi1 = max(yi);
ci = yi*(xi1-xi0)+xi;
[ci,num] = sort(ci);
xi = xi(num); yi = yi(num);
if floor(li/2)~=li/2
        xi = [xi xi(li)];
        yi = [yi yi(li)];
end

% Organize to pairs and separate by  NaN's
li = length(xi);
xi = reshape(xi,2,li/2);
yi = reshape(yi,2,li/2);
xi = [xi; ones(1,li/2)*nan];
yi = [yi; ones(1,li/2)*nan];
xi = xi(:)'; yi = yi(:)';

% Transform to the original coordinates
yi = yi*step;
xy = xi*ca-yi*sa;
yi = xi*sa+yi*ca;
xi = xy/xsc+x0;
yi = yi/ysc+y0;

hl = line('xdata',xi,'ydata',yi,...
    'color',color,'linestyle',style,'linewidth',width);
set(gca,'units',oldu)   % Set axes units back

if nargout>0, h=hl; end