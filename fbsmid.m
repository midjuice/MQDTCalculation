function [fmid,gmid] = fbsmid(r,element,channel,ll,B,gJ,Nmax,beta6)
    
    [bjpmat,bjmmat] = bmCal(element,channel,ll,B,Nmax,gJ);
    [~,v] = deltaCal(element,channel,ll,B,gJ);


    [fv,gv] = fbarinf(r,beta6,v,bjpmat,bjmmat);
    for i = 1:length(v)
        fmid(i) = cos(pi*v(i)/2).*fv(i) + sin(pi*v(i)/2).*gv(i);
        gmid(i) = -sin(pi*v(i)/2).*fv(i) + cos(pi*v(i)/2).*gv(i);
    end

    fmid = fmid(2:end);
    gmid = gmid(2:end);
end