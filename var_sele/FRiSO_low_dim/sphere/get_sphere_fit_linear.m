% @title Individually penalized ridge Frechet regression for Spherical data equipped with geodesic distance  
% @param Y: nx3 spherical output data
% @param xin: nxp potential Euclidean predictors
% @param xout: n2xp Euclidean predictors to predict sphere from

function [lc_fit,df] = get_sphere_fit_linear(Y, xin, xout)

[n,p]=size(xin);
  lc_fit = zeros(3, size(xout,1));
  M = spherefactory(3);
  ops.verbosity = 0;
  fr.M = M; lc.M = M;
  xinbar=mean(xin);
  % Covariance matrix weighted by lambda
  temp=inv(cov(xin)); %diag(sqrt(lam))*inv(diag(sqrt(lam))*cov(xin)*diag(sqrt(lam))+eye(p)*(1/n))*diag(sqrt(lam));
  % Weight matrix
  Kmat=1+(xin-repmat(xinbar, n,1))*temp*((xout-repmat(xinbar, size(xout,1),1))');
  
    if(nargout==2)
      if(max(max(abs(xin-xout)))>1e-6)
          error('Xout is not same as Xin; cannot output df!!!');
      end
      df=sum(diag(Kmat)./sum(Kmat,2));
    end

  for j = 1:length(xout)
    % initial guess (extrinsic) for trustregions algorithm
    y0 = sum(cell2mat(arrayfun(@(k) Kmat(k, j)*Y(:, k), 1:n, 'UniformOutput', false))')'; y0 = y0/norm(y0); 

    if(length(find(Kmat(:, j))) < 2)

      lc_fit(:, j) = repmat(NaN, 1, 3);

    else
      % Collect cost, gradient and Hession functions
      % Kmat is weights in front of the distance function, Y is the spherical data, 
      % M is the spherical object which defines what the optimizer should be and what the distance is
      lc.cost = @(y) get_cost(Kmat(:, j)', Y, y, M);
      lc.egrad = @(y) get_egrad(Kmat(:, j)', Y, y, M);
      lc.ehess = @(y, u) get_ehess(Kmat(:, j)', Y, y, M, u);

      % perform optimization 
      lc_fit(:, j) = trustregions(lc, y0, ops);

    end
    
  end

end
