% title: Modified coordinate descent algorithm for spherical data
% @param Y: training spherical data
% @param xin: training predictor data
% @param tau: specified tau
% @param M: spherical object for Manopt trustregions function

function lambdacur=sphere_fit_codesc(Y,xin,tau, M)

[n,p]=size(xin);

lambdacur=repmat(tau/p,1,p);
curobj=sphere_fit_codesc_obj(Y,xin,lambdacur, M);

wu=1;
mycount=1;
eps0=1e-4;
  
tttcand=[0 0.25 0.5 0.75 1];

% runs as long as fewer than 20 iterations and difference between iterations not too small
while(wu)
    objold=curobj;
    lambdaold=lambdacur;
    for j=1:p
      lambdacurMj=lambdacur;
      lambdacurMj(j)=0;
       if(sum(lambdacurMj)>eps0)
        lambdacurMj=lambdacurMj/sum(lambdacurMj);
        lambdaj=0*lambdacur;
        lambdaj(j)=1;
        objcand=[];
        for kk=1:length(tttcand)
            ttt=tttcand(kk);
            objcand(kk)=sphere_fit_codesc_obj(Y,xin,((lambdaj-lambdacurMj)*ttt+lambdacurMj)*tau ,M);
        end
        [dump,indobj]=min(objcand);
        tttmin=tttcand(indobj);
        

        lambdacur=(lambdaj*tttmin+lambdacurMj*(1-tttmin))*tau;
       end % if
%       lambdacur
    end % for j=1:p
    
    if(max(abs(lambdacur-lambdaold))< tau*0.01) %1e-4) 
        wu=0;
    end

    if(mycount>20)
      wu=0;
      display('Takes more than 20 loops to converge!!!')
    end
    mycount=mycount+1;
end % while (wu)

end


