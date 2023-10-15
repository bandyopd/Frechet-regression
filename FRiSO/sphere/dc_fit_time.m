% compute fitting time for sphere data after selecting variables using dc
function dc_fit_time()
%load('result_test3.mat')
addpath(genpath('Manopt_2.0'))

% set sample size, correlation level, number of predictors, and noise level
n=100;
r=0.5;
p=9;
tau_noise=0.2;

% define m(x), the regression curve
reg_curveWu = @(t) [sqrt(1-t(:,5).^2).*cos(pi*(t(:,1)+t(:,9))), sqrt(1-t(:,5).^2).*sin(pi*(t(:,1)+t(:,9))), t(:,5)];

result = struct();
nsim = 100;

n_cores = str2num(getenv('SLURM_NTASKS'));
pool = parpool('local', n_cores);

parfor i = 1:nsim % run simulations in parallel
    display(i)
    rng(i)
    %   execute generation of data
    [Y,X, ~]=generate_data(n, reg_curveWu, p, r, tau_noise);
    
    % use dcov
    tic
    dcov = compute_dcov(X, Y.');
    result(i).dcov_time = toc;
    [~, ind] = sort(dcov, 'descend');
    % fitting time using 3 variables
    sele_ind = sort(ind(1:3));
    tic
    [~,~] = get_sphere_fit_linear(Y, X(:, sele_ind), X(:, sele_ind));
    result(i).fit_time1 = toc;
    % fitting time using 4 variables
    sele_ind = sort(ind(1:4));
    tic
    [~,~] = get_sphere_fit_linear(Y, X(:, sele_ind), X(:, sele_ind));
    result(i).fit_time2 = toc;

end
save('test3_dc_fit_time.mat', 'result');

delete(pool)
end