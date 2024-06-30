
function Eb = ebfunc(ll,B,mf1k,alpha1k,mf2k,alpha2k,gJ,gI,Ehf,Ia,Ks,Kt,E6,sgn)
    E = linspace(-100,0,1000);
    for i = 1:length(E)
        a = real(KeffuncID(E(i),ll,B,mf1k,alpha1k,mf2k,alpha2k,gJ,gI,Ehf,Ia,Ks,Kt,E6,sgn));
        if abs(a - tan((2*ll+1)*pi/8)) < 1
            Eb = E(i);
            break;
        else
            Eb = 1;
        end

    end







end