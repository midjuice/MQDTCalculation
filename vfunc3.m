function y=vfunc3(Delta)
dc=1.58082626;ll=3;Nmax=10;
if(abs(abs(Delta)-dc)<1e-8)
    nu=2-1e-9;
elseif(abs(Delta)<dc)
    nu=fzero(@(x) Lambdal(x,Delta,ll,Nmax),[1.72  2-1e-11]);
elseif(abs(Delta)<18)
    nu=2+1i*fzero(@(x) Lambdal(x*1i,Delta,ll,Nmax),[1e-8,2]);
elseif(abs(Delta)<102)
    nu=2+1i*fzero(@(x) Lambdal(x*1i,Delta,ll,Nmax),[1.9,4]);
else
    nu=2+1i*fzero(@(x) Lambdal(x*1i,Delta,ll,2*Nmax),[3.95,10]);
end

y=nu;
end