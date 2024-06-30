function [Xl,Yl,Ml,theta,Gv]=XYMfunc(v,Delta,ll,Nmax)
vo=(2*ll+1)/4;
[bjp,bjm,Nb,Cvp,Cvm]=bjfunc(v,Delta,ll,Nmax);
Xl=1;Yl=0;
for mm=1:fix(Nb/2)
    Xl=Xl+(bjp(2*mm)+bjm(2*mm))*(-1)^mm;
    Yl=Yl+(bjp(2*mm-1)-bjm(2*mm-1))*(-1)^(mm-1);
end
if(mod(Nb,2)>0)
    mm=fix(Nb/2);
    Yl=Yl+(bjp(2*mm+1)-bjm(2*mm+1))*(-1)^mm;
end
atemp=gammacomplex(1-v)*gammacomplex(-v+vo+1)*gammacomplex(-v-vo+1);
btemp=gammacomplex(v+1)*gammacomplex(v-vo+1)*gammacomplex(v+vo+1);
Ml=(abs(Delta))^(2*v)*(Cvm/Cvp)*(atemp/btemp);
theta=atan(Yl/Xl);
if(Xl<0)
    theta=theta+pi;
end
Gv=(abs(Delta))^(-v)*gammacomplex(1+vo+v)*gammacomplex(1-vo+v)*Cvp/gammacomplex(1-v);
end