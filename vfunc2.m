function y=vfunc2(Delta)
dc=0.4306920854;ll=2;Nmax=20;
if(abs(abs(Delta)-dc)<1e-10)
    nu=1+1e-6;
elseif(abs(Delta)<dc)
    nu=fzero(@(x) Lambdal(x,Delta,ll,Nmax),[1+1e-8,1.25]);
elseif(abs(Delta)<14)
    nu=1+1i*fzero(@(x) Lambdal(x*1i,Delta,ll,Nmax),[1e-8,2]);
elseif(abs(Delta)<95)
    nu=1+1i*fzero(@(x) Lambdal(x*1i,Delta,ll,Nmax),[1.95, 4]);
else
    nu=1+1i*fzero(@(x)  Lambdal(x*1i,Delta,ll,2*Nmax),[3.98,10]);
end

y=nu;
end