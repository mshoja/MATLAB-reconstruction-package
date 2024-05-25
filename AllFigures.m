
% A short explanation for 3 data requirements

% There are three data requirements: 1) Data stationarity 2) Data
% Markovicity 3)Data resolution. If data are simulated then all we need to
% do is to check the third requirement since simulated data are stationary and
% Markov by their construction. However, if data are real then we have to check the
% three data requirements. 
% 
% If the first requirement fails, i.e, the dataset is not stationary then we should chunk data into smaller time periods
% so that the smaller datasets over these periods are stationary and repeat reconstruction over each time period. This
% is called 'moving-window' approach.
% 
% If the second requirement is not met, i.e., data are not Markov, then you should find a quntity
% called 'Markov-Einstein' (ME) time scale. If ME = 1 then your data is Markov and you can analyze the 
% entire dataset. However, if ME>1 then your data is not Markov but a rarified sample with every ME-th data points
% is Markov and analysis must be performed on this rarified sample (or more rarified samples) not the entire dataset.

% The third requirement (data resolution):
% For all datasets (real or simulated) it is necessay to check data resolution and find to
% which category data resolution falls. As explained for the second requirement, you must check the
% resolution of a rarified sample which is Markov not the entire dataset. To investigate data resolution 
% you should find a quantity called the 'Relaxation time' (RT) using the command 'RelaxationTime'. Here, we
% have 3 categories as below:

% 1) RT<=50, which signifies that data falls in the category of 'low-resolution'
% 1) 50<RT<=100, which signifies that data falls in the category of 'medium-resolution'
% 1) RT>100, which signifies that data falls in the category of 'high-resolution'
% If data resolution is medium or high then you can follow 'Euler'
% reconstruction. Otherwise a 'Hermite' reconstruction must be followed

% We explain these data requirements as well as how to reconstruct a dataset in details across the 5 datasets
% corresponding with Figures 1-5 in the paper. Two datasets are real datasets from Ecology (Figure 3) and the
% climate(Figure 5). The details of these datasets are discussed in the paper.

%%
% Figure 1 

% This is a simulated dataset. So, we only check the third data requirement in below
S = load('OUdata1D.mat');data = S.data;  %load the data
data=data(1:20000); %This is a big timeseries with 10^6 data points. Here, we just use 
%its first 20000 data points
RT = RelaxationTime(data);

% Here, we find that RT = 94.4225 (50<RT<100) meaning that this dataset has a medium 
% (but close to high) resolution. So, we follow a Euler reconstruction in below

%******** A) Fitting a parametric model to data
% Reconstructing the data
dt=0.01;   % this is the time step we used to generate this dataset. Note that it can be arbitrary but we used the same dt for the reprodicibility of the true parameters
mu = @(x,par)par(1).*x;sigma = @(x,par)par(2);   % since mu and sigma are function handles this means that we want to consider parametric modeling
result1 = euler_reconstruction(data, dt, 'mu', mu, 'sigma', sigma,  ...
'gradient_fun', eulergrad(mu, sigma), 'lb', [-200 eps], 'ub', [200 200]);

%******* B) Fitting a spline model to data
L=-2;R=1.8; %since we have a ‘spline’ model it is better to shrink the state space
data(data<L | data>R) = nan;  %This is VERY important: In spline modeling if you consider a smaller range for your data then you must assign ‘nan’ to those few data points falling outside this range.
mu=8;sigma=8;  %since mu and sigma are numbers this means that we want to consider spline modeling with 8 knots for mu and 8 knots for sigma
result2 = euler_reconstruction(data, dt, 'nKnots', [mu sigma], 'spline', 'CC', 'L', ...
L, 'R', R, 'lb', [zeros(1, mu) - 10, zeros(1, sigma)+eps], 'ub', zeros(1, mu + sigma) + 10, 'solver', 'fmincon'); % 'CC' means that we want to consider cubic splines for both drift and diffusion functions

%********Plotting Figure1
% Left panel
S = load('OUdata1D.mat');data = S.data;data=data(1:20000);
t = linspace(0,200,length(data));
plot(t,data,'-k');
xlabel('time [a.u]');
ylabel('state [a.u]');
title('Ornstein-Uhlenbeck model');

% Right panel
mu=@(x,par)par(1).*x;sigma=@(x,par)par(2)+0.*x; %this is true model  
par=zeros(2,1);par(1)=-1;par(2)=1; %true model parameters
xplot=linspace(L,R,1000); % a dense mesh across the considered range
plot_results(result2,xplot,mu(xplot,par),sigma(xplot,par));

%%

% Figure 2

% This is a simulated dataset. So, we only check the third data requirement in below
S = load('MayData1D.mat');data = S.data;  %load the data
data=data(1:100000); %This dataset has 3*10^5 data points. We only use the first third of the dataset 
RT = RelaxationTime(data);

% Here, we find that RT = 3.9785e+03 meaning that this dataset has a high resolution. So, we follow a Euler reconstruction

%******** A) Fitting a parametric model to data
% NOTE: % Since this model is nonlinear you should run the following code lines a few times (Alternatively, you can consider increasing the
% 'search_agents' from 5 (default) to a bigger value) to make sure you do not miss the global solution. There is a possibility to get different
% outcomes each time. But, the solution with the lowest objective function (- sum of log-likelihoods) is the fitest solution. These information apear on the command window  

dt=0.01;  % This is the time step we used tosimulate this dataset (but it is arbitrary)
mu = @(x,par)par(1).*x.*(1-x./par(2))-par(3).*x.^2./(x.^2+par(4).^2);sigma=@(x,par)par(5);
result1 = euler_reconstruction(data, dt, 'mu', mu, 'sigma', sigma, 'gradient_fun',  eulergrad(mu, sigma), ...
    'lb', zeros(1,5)+eps, 'ub', 15.*ones(1,5),'useparallel',true,'search_agents', 5);  

%******* B) Fitting a spline model to data
dt = 0.01;
L = min(data);R = max(data);  % here, the entire range of dataset is considered
mu = 8;sigma = 8; % A spline model with 8 knots for mu and 8 knots for sigma
result2 = euler_reconstruction(data, dt, 'nKnots', [mu sigma], 'spline', 'CC', 'L', ...
L, 'R', R, 'lb', [zeros(1, mu) - 10, zeros(1, sigma)+eps], 'ub', zeros(1, mu + sigma) + 10, 'solver', 'fmincon', 'search_agents', 1);   % 'CC' means that we want to consider cubic splines for both drift and diffusion functions

%********Plotting Figure2
% Left panel
plot(data,'-k');
xlabel('time');ylabel('state, x');
title('Grazing model of May');

% Right panel
r=1.01;K=10;g=2.75;a=1.6;s=0.4;  % true parameter values
par = [r K g a s]; 
mu = @(x,par)r.*x.*(1-x./K)-g.*x.^2./(x.^2+a.^2);sigma=@(x,par)s;  % true model
xplot=linspace(L,R,2000);
plot_results(result2,xplot,mu(xplot,par),sigma(xplot,par));

%%

% Figure 3

% This is a real dataset (phycocyanin levels in Lake Mendota). So, we should check the three data requirement in below
data = readmatrix('BGA_stdlevel_2011.csv');data = data(:,3);   % Load the data

%1) Test of stationarity (if the test result is 1, as it is, then data are stationary)
[~,~,~,~,reg] = adftest(data,'model','ARD','lags',0:15);  % the input 0:15 is the number of lags we try in fitting an autoregressive model to data (it should be big enough. See next command). 
[~,lagndx] = min([reg(:).BIC]);  % this tells us how many lags we need (lagndx is 8<15. So, considering 15 lags in the previous command was enough. Otherwise, increase it and repeat the process)
[h, Pvalue,~]=adftest(data,'model','ARD','lags',lagndx);   % h is the test result with a pvalue of Pvalue (here h = 1, so the dataset is stationary)

%2) Test of Markovicty 
%NOTE: Here you need the 'Burg' folder. So, you should add its path to your
%working folder by typing the command 'addpath('The path of the Burg folder')' 
order = 10;  % Number of autoregressive lags to consider (should be big enough. Continue reading ... you will understand if this lag (10) is enough or not)
AR = ME_TimeScale(data,order);

% Here you get the following AR vector
%1.0000   -0.6261   -0.2177   -0.0840   -0.0339   -0.0120   -0.0068   -0.0009   -0.0080    0.0005   -0.0089

% Ignore the first element (which is always 1). The lag after which the
% rest of lements of AR are 0 is the 'Markov-Einstein' (ME) time scale. This real dataset shows long-range correlations. However, 
% after lag 3 (-0.0840) the correlations are very small and we can roughly say that ME time scale is 3. A dataset is Markov if it has a ME time scale of 1. Otherwise, the
% entire dataset is not Markov. However, in this dataset if we consider every third data point, data(1:3:end), then this rarified sample is Markov. We call it Markov sample 

% Checking the resolution of Markov sample
% Note: We must investigate the resolution of a rarified sample which is Markov not the dataset. Here the Markov
% sample consist of every third data point data(1:3:end).
RT = RelaxationTime(data(1:3:end));   % Relaxation time of the Markov sample

% Here, we find that RT = 273.7843>100 meaning that the Markov sample has a high resolution. So, we apply a Euler reconstruction

% Now, we apply spline reconstruction to the Markov sample data(1:3:end)
data = data(1:3:end); % As explained, a rarified sample with every third data point is Markov
dt = 1;  % This is completely arbitrary. Note that in spline modeling, (since splines are linear functions of parameters) if you, for instance
% would consider dt = 2 then drift parameters will be divided by 2 and diffusion parameters will be divided by sqert(2), so on.  
L = -6.5;R = 6;   % Here we shrink the range of data a bit since there is very few data points near the data borders
mu = 8;sigma = 8; % A spline model with 8 knots for mu and 8 knots for sigma
result = euler_reconstruction(data, dt, 'nKnots', [mu sigma], 'spline', 'CC', 'L', ...
L, 'R', R, 'lb', [zeros(1, mu) - 10, zeros(1, sigma)+eps], 'ub', zeros(1, mu + sigma) + 10, 'solver', 'fmincon', 'search_agents', 5);

% Left panel
A = readmatrix('BGA_stdlevel_2011.csv');
T = 160:10:250;
T1 = zeros(1,length(T));
time = A(1:end-8,2);
for i=1:length(T)
    T1(i) = find(time>T(i),1);
end
my_color = [0 0.4470 0.7410];
plot(data,'color',my_color);
xticks(T1);
xticklabels(T);
xlim([-inf inf]);
ylim([-inf inf]);
xlabel('Time (day number in 2011)');
ylabel('Phycocyanin level');
title('Lake Mendota');

% Right panel
xplot=linspace(L,R,2000);
plot_results(result,xplot)

%%

% Figure 4

% This is a simulated dataset. So, we only check the third data requirement in below
S = load('MayData1D.mat');data = S.data;  %load the data
data=data(1:300:end);  %We consider every 300-th data point (this is to generate a dataset which is low-resolution and show that Euler reconstruction does not suffice and a Hermite reconstruction is needed)
RT = RelaxationTime(data);

% Here, we find that RT = 12.7959<50 meaning that this dataset has a low resolution. So, we must follow a Hermite reconstruction reconstruction

% Applying Hermite reconstruction using spline modeling 
% This has two phases: in the first phase we still follow Euler reconstruction to get a rough result. In the second phase a Hermite
% reconstruction is performed to improve the Euler reconstruction outcomes

% First and second phases (Euler and Hermite reconstructions)
S = load('MayData1D.mat');data = S.data;   
data=data(1:300:end); 
L = 0;R = max(data); %we consider the entire range of data
% data(data<L | data>R) = nan;  
dt = 3;  % Since in the mother dataset dt = 0.01 and here we considered every 300-th data points the actual time step is 300*0.01 = 3  
mu = 8; sigma = 8;  
% AN EXTREMELY IMPORTANT POINT: To improve execution time and also accuracy we strongly recommend to use 'quardatic' splines 
% using the flag 'QQ' when Hermite reconstruction is needed, as is the case here
result = euler_reconstruction(data, dt, 'nKnots', [mu sigma], 'spline', 'QQ', 'L', L, 'R', R, ..., 
'lb', [zeros(1, mu) - 10, zeros(1, sigma)+eps], 'ub', zeros(1, mu + sigma) + 10, 'solver', 'fmincon', 'search_agents', 5);
legpoints = legitimate_points(data, dt, 'prev', result, 'prev_range', 0.5, 'j', 3, 'k', 9);   % Note that for low-resolution data a K value between 6 and 12 is needed (here we used K = 9)
result_her = hermite_reconstruction(data, dt, 'prev', legpoints, 'solver',...,
    'fmincon');

% Left panel
mu = @(x,par)par(1).*x.*(1-x./par(2))-par(3).*x.^2./(x.^2+par(4).^2);sigma=@(x,par)par(5)+0.*x;  % true model
par = zeros(5,1);par(1) = 1.01;par(2)= 10;par(3) = 2.75;par(4) = 1.6;par(5) = 0.4; %true model parameters
xplot=linspace(L,R,2000); % a dense mesh across the considered range
plot_results(result,xplot,mu(xplot,par),sigma(xplot,par)); % Euler & true models 

% Right panel
plot_results(result_her,xplot,mu(xplot,par),sigma(xplot,par));

%%

% Figure 5

% This is a real dataset (a climate record covering the northern hemisphere climate from 70000 to 20000 years before the present time). So, we must check the three data requirement in below
data = readmatrix('NGRIP20.csv');   % Load the data
data=data(2649:5081);  % this corresponds with time period from 70000 to 20000 years before the present time

%1) Test of stationarity (if the test result is 1, as it is, then data are stationary)
[~,~,~,~,reg] = adftest(data,'model','ARD','lags',0:10);  % the input 0:10 is the number of lags we try in fitting an autoregressive model to data (it should be big enough. See next command). 
[~,lagndx] = min([reg(:).BIC]);  % this tells us how many lags we need (lagndx is 5<10. So, considering 10 lags in the previous command was enough. Otherwise, increase it and repeat the process)
[h, Pvalue,~]=adftest(data,'model','ARD','lags',lagndx);   % h is the test result with a pvalue of Pvalue (here h = 1, so the dataset is stationary)

%2) Test of Markovicty 
%NOTE: Here you need the 'Burg' folder. So, you should add its path to your
%working folder by typing the command 'addpath('The path of the Burg folder')' 
order = 10;  % Number of autoregressive lags to consider (should be big enough. Continue reading ... you will understand if this lag (10) is enough ro not)
AR = ME_TimeScale(data,order);

% Here you get the following AR vector
%1.0000   -0.4796   -0.2098   -0.1082   -0.0616   -0.0197   -0.0468

% Ignore the first element (which is always 1). The lag after which the
% rest of lements of AR are 0 is the 'Markov-Einstein' (ME) time scale. This real dataset shows long-range correlations. However, 
% after lag 3 (-0.1082) the correlations are rather samll and we can  roughly say that ME time scale is 3. Unfortunately, this dataset is quite
% small with only 2433 data points, so considering a rarified sample with every third points severely shrinks the data. We, therefore, decided to
% consider every othewr data points data(1:2:end) for the analysis, so the ME = 2.  A dataset is Markov if it has a ME time scale of 1. Otherwise, the
% entire dataset is not Markov. However, in this dataset if we consider every other data point, data(1:2:end), then this rarified sample is nearly Markov. We call it Markov sample 

% Checking the resolution of Markov sample
% Note: We must investigate the resolution of a rarified sample which is Markov not the dataset. Here the Markov
% sample consist of every third data point data(1:3:end).
RT = RelaxationTime(data(1:2:end));   % Relaxation time of the Markov sample

% Here, we find that RT = 15.0951<50 meaning that the Markov sample has a low-resolution. So, we apply a Hermite reconstruction

% Applying Hermite reconstruction using spline modeling 
% This has two phases: in the first phase we still follow Euler reconstruction to get a rough result. In the second phase a Hermite
% reconstruction is performed to improve the Euler reconstruction outcomes

% First and second phases (Euler and Hermite reconstructions)
data = readmatrix('NGRIP20.csv');
data=data(2649:5081);
data = data(1:2:end);  % Markov sample
L = -45.5;R = -38.2;  % We shrink the data range a bit since near data borders there is very few data points (which adversely impact the analysis)
dt = 1;  % This is completely arbitrary. Note that in spline modeling, (since splines are linear functions of parameters) if you, for instance
% would consider dt = 2 then drift parameters will be divided by 2 and diffusion parameters will be divided by sqert(2), so on.
data(data<L | data>R) = nan;
mu = 7; sigma = 7;  
% AN EXTREMELY IMPORTANT POINT: To improve execution time and also accuracy we strongly recommend to use 'quardatic' splines 
% using the flag 'QQ' when Hermite reconstruction is needed, as is the case here
result = euler_reconstruction(data, dt, 'nKnots', [mu sigma], 'spline', 'QQ', 'L', L, 'R', R, ..., 
'lb', [zeros(1, mu) - 10, zeros(1, sigma)+eps], 'ub', zeros(1, mu + sigma) + 10, 'solver', 'fmincon', 'search_agents', 20); % we used 20 search agents (default is 5)
legpoints = legitimate_points(data, dt, 'prev', result, 'prev_range', 0.5, 'j', 3, 'k', 9);   % Note that for low-resolution data a K value between 6 and 12 is needed (here we used K = 9)
result_her = hermite_reconstruction(data, dt, 'prev', legpoints, 'solver', 'fmincon');

% Top panel
data = readmatrix('NGRIP20.csv');
data=data(2649:5081);
data = data(1:2:end);
my_color = [0 0.4470 0.7410];
plot(2649:2:5081,data,'-','color',my_color);
age = linspace(2649,5081,6);
xticks(age);xticklabels([70 60 50 40 30 20]);
yticks(-47:2:-37);yticklabels(-47:2:-37);
ylabel('\delta^{18}O (Permille)');
xlabel('Age (ky before present)');
title('A low-resolution climate recrord');
xlim([-inf inf]);ylim([-47 -36]);

% bottom panels
xplot=linspace(L,R,2000); % a dense mesh across the considered range
plot_results(result,xplot);
plot_results(result_her,xplot,'eff_potential');  % Here, we can also ask the packge to plot the effective potential (it is similar to the usual concept of potential function).
