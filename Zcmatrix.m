function [Zfs,Zfc,Zgs,Zgc]=Zcmatrix(Es,ll)
Nmax=10;
Delta=Es/16;
v=vfunctot(Delta,ll);
vo=(2*ll+1)/4;
[Xl,Yl,Ml,~,Gv]=XYMfunc(v,Delta,ll,Nmax);
Al=Gv*cos(pi*(v-vo))/(sqrt(2)*(Xl^2+Yl^2)*sin(pi*v));

tanalpha=tan(pi*(v-vo));sinv=sin(pi*v/2);cosv=cos(pi*v/2);
Zfs=((1-(-1)^ll*Ml*tanalpha)*sinv*Xl+(1+(-1)^ll*Ml*tanalpha)*cosv*Yl)*Al;
Zfc=((tanalpha-(-1)^ll*Ml)*sinv*Xl+(tanalpha+(-1)^ll*Ml)*cosv*Yl)*Al;
Zgs=(-(1-(-1)^ll*Ml*tanalpha)*sinv*Yl+(1+(-1)^ll*Ml*tanalpha)*cosv*Xl)*Al;
Zgc=(-(tanalpha-(-1)^ll*Ml)*sinv*Yl+(tanalpha+(-1)^ll*Ml)*cosv*Xl)*Al;

%detZ=det([Zfs,Zfc;Zgs,Zgc])-1

end