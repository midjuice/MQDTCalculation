function y=vfunc1(Delta)
dc=0.147379179523;ll=1;Nmax=10;
if(abs(abs(Delta)-dc)<1e-10)
    nu=1-1e-6;
elseif(abs(Delta)<dc)
    nu=fzero(@(x) Lambdal(x,Delta,ll,Nmax),[0.75,1-1e-8]);
elseif(abs(Delta)<12)
    nu=1+1i*fzero(@(x) Lambdal(x*1i,Delta,ll,Nmax),[1e-8,2]);
elseif(abs(Delta)<91)
    nu=1+1i*fzero(@(x) Lambdal(x*1i,Delta,ll,Nmax),[1.9,4]);
else
    nu=1+1i*fzero(@(x) Lambdal(x*1i,Delta,ll,2*Nmax),[3.98,10]);
end

y=nu;
end