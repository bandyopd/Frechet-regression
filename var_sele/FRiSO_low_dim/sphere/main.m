% @title Main variable selection for spherical data
% @description This function performs FRiSO for the spherical data simulation in example 5.4.1
% @details:
% 1. Sets data generation details
% 2. Defines m(x), the regression curve
% 3. Generates validation data and testing data
% 4. Plots the generated data
% 5. 100 times, generates training data, given tau, performs modified coordinate descent, 
%    selects the best tau based on validation error, then collects testing error based on best model

addpath(genpath('Manopt_2.0'))

% set sample size, correlation level, number of predictors, and noise level
n=200;
r=0.5;
p=9;
tau_noise=0.2;

% define m(x), the regression curve
reg_curveWu = @(t) [sqrt(1-t(:,5).^2).*cos(pi*(t(:,1)+t(:,9))), sqrt(1-t(:,5).^2).*sin(pi*(t(:,1)+t(:,9))), t(:,5)];


% generate validation data for picking optimal tau
rng(2022)
[Ytune,Xtune, ~]=generate_data(n, reg_curveWu, p, r, tau_noise);

result = struct();
nsim = 100;
taus=(1:2:80)*0.1*0.5;

% Set up parallel processing

%pool = parpool(2);
%  display('Where am I 6');
n_cores = str2num(getenv('SLURM_NTASKS'));
pool = parpool('local', n_cores);

parfor i = 1:nsim % run simulations in parallel
    rng(i)
    %   execute generation of data
    [Y,X, ~]=generate_data(n, reg_curveWu, p, r, tau_noise);
    
    % use dcov
    dcov = compute_dcov(X, Y.');
    % 100 times, for each tau, perform coordinate descent, validation error, best tau, and testing error
    [lambdacur,RSS,df,tuneer]=forjloop(Y,X,taus, Ytune, Xtune);
    result(i).X=X;
    result(i).Y=Y;
    result(i).lambdacur=lambdacur;
    result(i).RSS=RSS;
    result(i).df=df;
    result(i).tuneer=tuneer;
    result(i).dcov = dcov;
end
save('resultLowNoise.mat', 'result');

delete(pool)

