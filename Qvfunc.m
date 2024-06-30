function y=Qvfunc(v,Delta,ll,Nmax)
vo=(2*ll+1)/4;
fj=@(jj) 1/((v+jj)*((v+jj)^2-vo^2)*(v+jj+1)*((v+jj+1)^2-vo^2));
Qv=1;
for kk=0:Nmax-1
    Qv=1/(1-Delta^2*fj(Nmax-kk)*Qv);
end
y=Qv;
end