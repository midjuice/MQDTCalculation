function [cj,Nc]=cjfunc(v,Delta,ll,Nmax)
qvj=Qvfunc(v,Delta,ll,Nmax);
cj=[0 0];
cj(1)=qvj;Nc=1;
while (abs(1-qvj)>1e-8)
    Nc=Nc+1;
    qvj=Qvfunc(v+Nc-1,Delta,ll,Nmax);
    cj(Nc)=qvj*cj(Nc-1);
end

end