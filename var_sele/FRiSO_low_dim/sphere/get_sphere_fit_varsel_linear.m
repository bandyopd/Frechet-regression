% @title: Individually penalized ridge Frechet regression for spherical data
% @descriptions: Gets sphere fits for simulation data
% @param Y: training spherical data
% @param xin: training predictor data
% @param xout: data to be predicted on
% @param lam: given individual penalty weights, lambda 
% @details: 
% 1. Prepares data for predictions
% 2. Computes weighted inverse matrix
% 3. Computes cost, gradient and Hessian functions given data
%    Manopt defaults to the geodesic distance when the data is spherical
% 4. Performs trustregions optimization from Manopt for fitted values

function [lc_fit,df] = get_sphere_fit_varsel_linear(Y, xin, xout, lam)

[n,p]=size(xin);
  lc_fit = zeros(3, size(xout,1));
  M = spherefactory(3);
  ops.verbosity = 0;
  fr.M = M; lc.M = M;
 
  xinbar=mean(xin);

  % Computes weighted inverse matrix within the cost function
  temp=diag(sqrt(lam))*inv(diag(sqrt(lam))*cov(xin)*diag(sqrt(lam))+eye(p)*(1/n))*diag(sqrt(lam));
  Kmat=1+(xin-repmat(xinbar, n,1))*temp*((xout-repmat(xinbar, size(xout,1),1))');
  
    if(nargout==2)
      if(max(max(abs(xin-xout)))>1e-6)
          error('Xout is not same as Xin; cannot output df!!!');
      end
      df=sum(diag(Kmat)./sum(Kmat,2));
    end
  
  for j = 1:size(xout,1) %length(xout)
    y0 = sum(cell2mat(arrayfun(@(k) Kmat(k, j)*Y(:, k), 1:n, 'UniformOutput', false))')'; y0 = y0/norm(y0); % Initial guess (extrinsic)
    if(length(find(Kmat(:, j))) < 2)
      lc_fit(:, j) = repmat(NaN, 1, 3);

    else
      % computes cost, gradient and Hessian function for trustregions function
      lc.cost = @(y) get_cost(Kmat(:, j)', Y, y, M);
      lc.egrad = @(y) get_egrad(Kmat(:, j)', Y, y, M);
      lc.ehess = @(y, u) get_ehess(Kmat(:, j)', Y, y, M, u);
      % gets fitted values from nonlinear optimizer, trustregions (Manopt)
      lc_fit(:, j) = trustregions(lc, y0, ops);

    end
    
  end

end
