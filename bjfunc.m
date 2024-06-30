function [bjp,bjm,Nb,Cvp,Cvm]=bjfunc(v,Delta,ll,Nmax)
vo=(2*ll+1)/4;
[cjp,Ncp]=cjfunc(v,Delta,ll,Nmax);
[cjm,Ncm]=cjfunc(-v,Delta,ll,Nmax);
dmt=gammacomplex(v)*gammacomplex(v-vo+1)*gammacomplex(v+vo+1);
dnt=gammacomplex(v+1)*gammacomplex(v-vo)*gammacomplex(v+vo);
fp=@(jj) (-Delta)^jj*dmt/(gammacomplex(v+jj)*gammacomplex(v-vo+1+jj)*gammacomplex(v+vo+1+jj));
fm=@(jj) (-Delta)^jj*gammacomplex(v-jj+1)*gammacomplex(v-vo-jj)*gammacomplex(v+vo-jj)/dnt;
fjp=fp(1);fjm=fm(1);
bjp=[0 0];bjm=[0 0];
bjp(1)=fjp*cjp(1);bjm(1)=fjm*cjm(1);
Nb=1;
while (abs(fjp)+abs(fjm)>1e-17)
    Nb=Nb+1;
    if(Nb<Ncp+1)
        cjpjj=cjp(Nb);
    else
        cjpjj=cjp(Ncp);
    end
    
    if(Nb<Ncm+1)
        cjmjj=cjm(Nb);
    else
        cjmjj=cjm(Ncm);
    end
    fjp=fp(Nb);fjm=fm(Nb);
    bjp(Nb)=fjp*cjpjj;bjm(Nb)=fjm*cjmjj;
    
end
Cvp=cjp(Ncp);Cvm=cjm(Ncm);
end
