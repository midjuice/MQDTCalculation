function y=vfunctot(Delta,ll)
if(ll==0)
    y=vfunc0(Delta);
elseif(ll==1)
    y=vfunc1(Delta);
elseif(ll==2)
    y=vfunc2(Delta);
elseif(ll==3)
    y=vfunc3(Delta);
elseif(ll==4)
    y=vfunc4(Delta);
elseif(ll==5)
    y=vfunc5(Delta);
else
    y=NaN;
    warning='ll is out of range!'
    return
end

end