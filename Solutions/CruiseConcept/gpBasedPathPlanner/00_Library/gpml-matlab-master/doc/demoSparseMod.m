disp('See http://www.gaussianprocess.org/gpml/code/matlab/doc/ for details.')
clear all, close all, write_fig = 0; N = 20;
sd = 3; rand('seed',sd), randn('seed',sd)       % set a seed for reproducability

fprintf('a) switch between FITC/VFE/SPEP via the opt.s parameter\n')
a = 0.3; b = 1.2; f = @(x) a*x + b + sin(x);               % underlying function
n = 30; sn = 0.5;          % number of training points, noise standard deviation
x = 2*rand(n,1)-1; 
x = 1+4*x+sign(x); 
y = f(x)+sn*randn(n,1);      % sample data

cov = {@covSEiso}; sf = 2; ell = 1.0; hyp.cov = log([ell;sf]);
mean = {@meanSum,{@meanLinear,@meanConst}}; 
hyp.mean = [a;b];
lik = {@likGauss};    
hyp.lik = log(sn); 
inf = @infGaussLik;

fprintf('Optimise hyperparameters.\n')
hyp = minimize(hyp,@gp,-N,inf,mean,cov,lik,x,y);      % optimise hyperparameters
xs = linspace(-8,10,2e3)'; ys = f(xs);                   % exact function values

[ymu,ys2] = gp(hyp,inf,mean,cov,lik,x,y,xs);                  % dense prediction

% turn into a sparse approximation
xu = linspace(-8,10,20)'; cov = {'apxSparse',cov,xu};           % inducing points
infv  = @(varargin) inf(varargin{:},struct('s',0.0));           % VFE, opt.s = 0
hyp.xu = xu;
hypInduced = minimize(hyp,@gp,-N,inf,mean,cov,lik,x,y);                % optimize with induced inputs as hyperparameter

[yu, yus] = gp(hypInduced,infv,mean,cov,lik,x,y, hypInduced.xu);                 % using induced inputs, generate induced outputs

cov = {@covSEiso}; % going back to original configuration, generate output with induced gp
[yInduced,sInduced] = gp(hypInduced,@infGaussLik,mean,cov,lik,hypInduced.xu,yu,xs);

confidenceBounds = [ymu+2*sqrt(ys2); flip(ymu-2*sqrt(ys2),1)];
confidencePoints = xs;
fill([confidencePoints; flip(confidencePoints)], confidenceBounds, 'y', 'DisplayName', '95% confidence - original',  'FaceAlpha', 0.2);
hold on;

confidenceBounds = [yInduced+2*sqrt(sInduced); flip(yInduced-2*sqrt(sInduced),1)];
confidencePoints = xs;
fill([confidencePoints; flip(confidencePoints)], confidenceBounds, 'g', 'DisplayName', '95% confidence - induced', 'FaceAlpha', 0.2);

plot(xs,ymu,'b-','LineWidth',2, 'DisplayName', 'original estimate'), hold on
plot(xs,ys,'k.','LineWidth',2, 'DisplayName', 'original data')
plot(xs,yInduced,'m:','LineWidth',2, 'DisplayName', 'induced based estimate')
plot(hypInduced.xu, yu,'r+', 'DisplayName', 'induction points');

xlim([-8,10]), ylim([-3,6]);
legend;
