function y=Lambdal(v,Delta,ll,Nmax)
vo=(2*ll+1)/4;
fp=@(v) Qvfunc(v,Delta,ll,Nmax)/((v+1)*((v+1)^2-vo^2));
y=v^2-vo^2-Delta^2*(fp(v)-fp(-v))/v;
end