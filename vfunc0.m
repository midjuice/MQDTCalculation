function y=vfunc0(Delta)
dc=0.09654417566;ll=0;Nmax=10;

if(abs(abs(Delta)-dc)<1e-10)
    nu=1e-6;
elseif(abs(Delta)<dc)
    nu=fzero(@(x) Lambdal(x,Delta,ll,Nmax),[1e-8,0.25]);
elseif(abs(Delta)<11)
    nu=1i*fzero(@(x) Lambdal(x*1i,Delta,ll,Nmax),[1e-8,2]);
elseif(abs(Delta)<88)
    nu=1i*fzero(@(x) Lambdal(x*1i,Delta,ll,Nmax),[1.97,4]);
else
    nu=1i*fzero(@(x) Lambdal(x*1i,Delta,ll,2*Nmax),[3.98,20]);
end

y=nu;
end