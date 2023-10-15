% @title Data generation for unit sphere data in R3 equipped with the geodesic distance

function [Y,X, mreg]=generate_data(n, reg_curve, p, r, tau)

% include correlation of predictors

halfcorr=chol(r.^abs(bsxfun(@minus,1:p,(1:p)')));
temp=randn(n,p)*(halfcorr);
X=normcdf(temp);
mreg=reg_curve(X)';
u=tau*randn(n,2);
Y = cell2mat(arrayfun(@(j) add_noise(mreg(:, j), u(j,:)), 1:n, 'UniformOutput', false));

end

