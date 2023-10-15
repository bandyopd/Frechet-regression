% @title: Objective function of cordinate descent algorithm
% @param Yin: training spherical data
% @param xin: training predictor data
% @param lam: current lambda value 
% @param M: spherical object specified for Manopt trustregions function
% @param xout: data to predict on

function [obj,lc_fit,df]=sphere_fit_codesc_obj(Yin,xin,lam, M, Yout, xout)

if(nargin==4)
    xout=xin;
    Yout=Yin;
end

if(nargout==3)
    [lc_fit,df] = get_sphere_fit_varsel_linear(Yin, xin, xout, lam);
else
    [lc_fit] = get_sphere_fit_varsel_linear(Yin, xin, xout, lam);
end

obj=sum(arrayfun(@(c) M.dist(Yout(:, c), lc_fit(:, c))^2, 1:size(Yout,2)));


end
