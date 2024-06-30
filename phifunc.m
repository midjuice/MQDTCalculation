function phi = phifunc(r,element,channel,ll,B,gJ,Nmax,beta6)
    Kc=Kcfunc(element,channel,B,ll,gJ);
    [f,g] = scatterfunc(r,element,channel,ll,B,gJ,Nmax,beta6);
    [h,j] = boundfunc(r,element,channel,ll,B,gJ,Nmax,beta6);
    ffinal = transpose([f,h]);
    gfinal = transpose([g,j]);
    phi = ffinal - Kc * gfinal;
    if B < 1040
        p = transpose([4.02486380827236e+84,2.74566739469417e+85,9.42336911798516e+85,2.71161639105832e+85,1.36049854403238e+86]);
    else
        p = transpose([1.59207486322476e+85,1.10080124806683e+86,4.24398692258415e+86,1.08748197610431e+86,6.04256201822567e+86]);
    end
    phi = phi./p;
end