function y=vfunc5(Delta)
dc=5.61577;ll=5;Nmax=10;
if(abs(abs(Delta)-dc)<1e-5)
    nu=2-1e-9;
elseif(abs(Delta)<dc)
    nu=fzero(@(x) Lambdal(x,Delta,ll,Nmax),[2.65 3-1e-11]);
elseif(abs(Delta)<29)
    nu=3+1i*fzero(@(x) Lambdal(x*1i,Delta,ll,Nmax),[1e-8,2]);
elseif(abs(Delta)<123)
    nu=3+1i*fzero(@(x) Lambdal(x*1i,Delta,ll,Nmax),[1.95,4]);
else
    nu=3+1i*fzero(@(x) Lambdal(x*1i,Delta,ll,2*Nmax),[3.95,10]);
end

y=nu;
end