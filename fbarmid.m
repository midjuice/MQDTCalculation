function [fv,gv] = fbarmid(r,beta6,v,bmp,bmm)
    fv = zeros(length(v),1);
    gv = zeros(length(v),1);
    for ii  = 1:length(v)
        [fv(ii,1),gv(ii,1)] =   fbarmidv(r,beta6,v(ii),bmp(ii,:),bmm(ii,:));
    end
end




function [ftot,gtot] = fbarmidv(r,beta6,v,bmp,bmm)
    fbtotp = 0;
    fbtotm = 0;
    gbtotp = 0;
    gbtotm = 0;
    for m = 1:length(bmp)
        [a,b] = fbar(r,beta6,v,bmp(m),m);
        [c,d] = fbar(r,beta6,v,bmm(m),-m);
        fbtotp = fbtotp + a;
        fbtotm = fbtotm + c;
        gbtotp = gbtotp + b;
        gbtotm = gbtotm + d;
    end
    ftot = fbtotp + fbtotm;
    gtot = gbtotm + gbtotp;
end






function [f,g] = fbar(r,beta6,v,bm,m)
    f = bm * r^(-1/2) * cbesselj(v+m,0.5*(r/beta6)^(-2));
    g = bm * r^(-1/2) * cbessely(v+m,0.5*(r/beta6)^(-2));
end