function y=vfunc4(Delta)
dc=2.073295;ll=4;Nmax=10;
if(abs(abs(Delta)-dc)<1e-6)
    nu=2-1e-9;
elseif(abs(Delta)<dc)
    nu=fzero(@(x) Lambdal(x,Delta,ll,Nmax),[2+1e-11 2.25]);
elseif(abs(Delta)<22)
    nu=2+1i*fzero(@(x) Lambdal(x*1i,Delta,ll,Nmax),[1e-8,2]);
elseif(abs(Delta)<111)
    nu=2+1i*fzero(@(x) Lambdal(x*1i,Delta,ll,Nmax),[1.9,4]);
else
    nu=2+1i*fzero(@(x) Lambdal(x*1i,Delta,ll,2*Nmax),[3.95,10]);
end

y=nu;
end