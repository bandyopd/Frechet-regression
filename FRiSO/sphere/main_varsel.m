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
  display('add path Manopt');

% set sample size, correlation level, number of predictors, and noise level
n=100;
r=0.5;
p=8;
tau_noise=0.2;

% define m(x), the regression curve
reg_curveWu = @(t) [sqrt(1-t(:,3).^2).*cos(pi*(t(:,5)+t(:,7))), sqrt(1-t(:,3).^2).*sin(pi*(t(:,5)+t(:,7))), t(:,3)];
  display('Where am I 3');

% generate validation data for picking optimal tau
rng(2022)
[Ytune,Xtune, ~]=generate_data(n, reg_curveWu, p, r, tau_noise);

% generate large testing data
[YtuneBig,XtuneBig, mregtuneBig]=generate_data(n*10, reg_curveWu, p, r, tau_noise);

% First plot the generated data
if 0   
    [s1 s2 s3] = sphere(128);
    lightGrey = 0.9*ones(1, 3);
    figure;
    set(gcf, 'PaperUnits', 'inches', 'PaperSize', [2.4, 2.6], 'PaperPosition', [0 0 2.4 2.6])
    surf(s1, s2, s3, 'FaceColor', lightGrey, 'EdgeColor', 'none', 'FaceAlpha', 0.4)
    view([1 1 1])
    hold on;
    plot3(mregtuneBig(1,:), mregtuneBig(2,:), mregtuneBig(3,:), 'kx', 'LineWidth', 0.5)
    plot3(YtuneBig(1,:), YtuneBig(2,:), YtuneBig(3,:), 'rx', 'LineWidth', 0.5)
end
  display('Where am I 4');

result = struct();
nsim=1
  display('Where am I 5');

% Set up parallel processing

%pool = parpool(2);
%  display('Where am I 6');

for i = 1:nsim % run simulations in parallel
    
%    clear lambdacur RSS df tuneer tuneBiger;
  display('Here 1');
    rng(i)
  display('Here 2');
%   execute generation of data
    [Y,X, ~]=generate_data(n, reg_curveWu, p, r, tau_noise);
  display('Here 3');
    taus=(1:2:80)*0.1*0.5;
  display('Here 4');

% 100 times, for each tau, perform coordinate descent, validation error, best tau, and testing error
[lambdacur,RSS,df,tuneer,tuneBiger]=forjloop(Y,X,taus, Ytune, Xtune, YtuneBig, XtuneBig);
  display('Here 6');  
  result(i).X=X;
  display('Here 7');
  result(i).Y=Y;
  result(i).lambdacur=lambdacur;
  result(i).RSS=RSS;
  result(i).df=df;
  result(i).tuneer=tuneer;
  result(i).tuneBiger=tuneBiger;
end
  save('resultLowNoise.mat');

%delete(pool)

