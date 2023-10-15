% @title Iterative function for FRiSO on Spherical data with the geodesic distance
% @description This function implements the coordinate descent algorithm given tau
% @param Y: an nx3 matrix of training spherical data outputs
% @param X: an nxp matrix of potential training data predictors from Euclidean space
% @param taus: a vector taus such that sum(lambdas)=tau
% @param Ytune: an n2 x 3 matrix of tuning spherical data outputs
% @param Xtune: an n2 x p matrix of tuning predictors
% @param YtuneBig: an n3 x 3 matrix of testing spherical data outputs
% @param XtuneBig: an n3 x 3 matrix of testing predictors
% @details:
% Loops through tau and perform optimization of lambdas utilizing modified coordinate
%    descent alogorithm

function [lambdacur,RSS,df,tuneer]=forjloop(Y,X,taus, Ytune, Xtune)
%addpath(genpath('/home/yichaowu/MatlabRun/NWshpereSimu1New/Manopt_2.0'))

% Define output is a sphere with geodesic distance
M = spherefactory(3);
% for each tau, optimize lambdas utilizing modified coordinate descent

for j=1:length(taus)
    %tic
    display(j)
    lambdacur(j,:)=sphere_fit_codesc(Y,X,taus(j), M);
    [RSS(j), ~, df(j)]=sphere_fit_codesc_obj(Y,X,lambdacur(j,:), M);
    tuneer(j)=sphere_fit_codesc_obj(Y,X,lambdacur(j,:), M, Ytune, Xtune);
    %toc
end
end

