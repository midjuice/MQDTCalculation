function [f,g] = boundfunc(r,element,channel,ll,B,gJ,Nmax,beta6)
    if r < 11.67
        [f,g] = fbsori(r,element,channel,gJ,ll,B,beta6);
    end
    if r > 11.67 && r < 11.735
        [f,~] = fbsmid(r,element,channel,ll,B,gJ,Nmax,beta6);
        [~,g] = fbsori(r,element,channel,gJ,ll,B,beta6);
    end
    if r > 11.735 && r < 655
        [f,g] = fbsmid(r,element,channel,ll,B,gJ,Nmax,beta6);
    end
    if r > 655
        [f,g] = fbsinf(r,element,channel,ll,B,gJ,Nmax,beta6);
    end
end