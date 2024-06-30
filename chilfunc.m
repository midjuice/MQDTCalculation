function y=chilfunc(Es,ll)
% if(Es>0)
%     msg = 'Error occurred. Input energy is greater than 0.';
%     error(msg)
% end

Delta=Es/16;
v=vfunctot(Delta,ll);
Nmax=10;
[Xl,Yl,Ml,~,~]=XYMfunc(v,Delta,ll,Nmax);
tantheta=Yl/Xl;
temp=tan(pi*v/2)*(1+Ml)/(1-Ml);
y=((tantheta+temp)/(1-tantheta*temp));

end