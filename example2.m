%model:
%mu=r*x*(1-x/K)*(x-C)
mu = @(x,par)par(1).*x.*(1-x./par(2)).*(x-par(3));
%sigma = s.*x
sigma = @(x,par)par(4).*x;

par = zeros(4, 1);
par(1) = 1; %r
par(2) = 10; %K
par(3) = 3; %C;
par(4) = 0.4; %s

x0 = 5;
deltat = 0.01;
L = 0.1;
R = 20;

%stabilize
y = simulate('parametric', L, R, @(x)mu(x,par), @(x)sigma(x,par), deltat, x0, 1000);

x0 = y(end);
nsteps=500000;
data = simulate('parametric', L, R, @(x)mu(x,par), @(x)sigma(x,par), deltat, x0, nsteps);

figure
ts=(1:nsteps).*deltat;
data=data(1:5:end);
ts=ts(1:5:end);

fprintf('The fitted relaxation time = %g\n', RelaxationTime(data))

plot(ts, data);
res = euler_reconstruction(data, deltat, 'mu', mu, 'sigma', sigma, 'ub', [10 20 10 10], 'lb', [0 0 0 0], 'search_agents', 5, 'solver', 'fmincon');
%res=euler_reconstruction(data,dt,'nknots',[4 1],'spline','CC');
xplot = linspace(min(data), max(data), 1000);
plot_results(res, xplot, mu(xplot, par), sigma(xplot, par));

% addpath 'C:\d\alg\MAT LAB\stats\LangvinToolbox\code'
%  % dt = 0.01;
%   M = 4; % M is the number of knots (do not choose too many knots). Note that for additive models the number of parameters is M+1=9. Too many knots mean too many parameters which should be estimated plus that at the end the results would be wobbly  
%   knots = linspace(min(data), max(data), M);
%   lb = [-10 .* ones(1, M) eps.* ones(1, 1)];
%   ub = [10 .* ones(1, M) 10 .* ones(1, 1)];
%  mesh_fine = linspace(knots(1), knots(end), 2000); %We need a fine mesh for our data range but please do not choose an extremely fine mesh. Here, the data range is [0 7.1678] (this is more than enough)
%  solver = 'fmincon';
%  ModelType = 'Additive noise';
%  %'L','Q','C','P','SCS'
%  SplineType = 'C'; %This means that we wish to fit a cubic spline as drift function
%  UseParallel = false;
%  Max_iter = realmax; %Unlike parametric models, for spline models we recommend to use ‘realmax’ as iteration number 
%  [estimated_par, f_best, mufun, sigmafun] = DiffusionReconst_Spline_Euler(data, dt, knots, mesh_fine, lb, ub, solver, ModelType, SplineType, [], Max_iter, UseParallel, 'Yes');
%  
