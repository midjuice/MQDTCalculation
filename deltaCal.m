function [Delta,v] = deltaCal(element,channel,ll,B,gJ)

    [Ia, mf1k, alpha1k, mf2k, alpha2k, gI, ~, Ehf, E6, ~, ~, ~] = Get_Params(element, channel);

    Es = 10e-8;
    Delta = zeros(1,1);
    v = zeros(1,1);
    Eth0 = (Ethfunc(B,mf1k(1),alpha1k(1),gJ,gI,Ia,Ehf)+Ethfunc(B,mf2k(1),alpha2k(1),gJ,gI,Ia,Ehf))/E6;
    for i = 1:length(mf1k)
        Ethk =  Ethfunc(B,mf1k(i),alpha1k(i),gJ,gI,Ia,Ehf)/E6+ Ethfunc(B,mf2k(i),alpha2k(i),gJ,gI,Ia,Ehf)/E6;
        Delta(i) = Es + Eth0 - Ethk;
        v(i)=vfunctot(Delta(i),ll);
    end
 
end