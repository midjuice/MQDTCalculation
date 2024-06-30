function det_value = bound_spectrum_eqn(E, ll, sgn, B, mf1k, alpha1k, mf2k, alpha2k, gJ, gI, Ia, Ehf, E6, Ks, Kt)
    Kc=zeros(1,1);
    x=(gJ-gI)*B*1.399624504/Ehf;
    Nc=length(mf1k);
    for kk=1:Nc
        for jj=1:Nc
            temp1 = Keffc(mf1k(kk),alpha1k(kk),mf2k(kk),alpha2k(kk),mf1k(jj),alpha1k(jj),mf2k(jj),alpha2k(jj),Ia,x,Ia,x,Ks,Kt);
            temp2 = Keffc(mf2k(kk),alpha2k(kk),mf1k(kk),alpha1k(kk),mf1k(jj),alpha1k(jj),mf2k(jj),alpha2k(jj),Ia,x,Ia,x,Ks,Kt);
            temp3 = Keffc(mf1k(kk),alpha1k(kk),mf2k(kk),alpha2k(kk),mf2k(jj),alpha2k(jj),mf1k(jj),alpha1k(jj),Ia,x,Ia,x,Ks,Kt);
            temp4 = Keffc(mf2k(kk),alpha2k(kk),mf1k(kk),alpha1k(kk),mf2k(jj),alpha2k(jj),mf1k(jj),alpha1k(jj),Ia,x,Ia,x,Ks,Kt);
            Kc(kk,jj) = (temp1+temp4+sgn*(-1)^ll*(temp3+temp2))/2;
            if (mf1k(kk)==mf2k(kk) && alpha1k(kk)==alpha2k(kk))
                Kc(kk,jj) = Kc(kk,jj)/sqrt(2);
            end
            if (mf1k(jj)==mf2k(jj) && alpha1k(jj)==alpha2k(jj))
                Kc(kk,jj) = Kc(kk,jj)/sqrt(2);
            end
        end
    end
    chi_array = zeros(1,Nc);
    Eth0 = Ethfunc(B,mf1k(1),alpha1k(1),gJ,gI,Ia,Ehf)/E6 + Ethfunc(B,mf2k(1),alpha2k(1),gJ,gI,Ia,Ehf)/E6;
    for i=1:Nc
        Eth = Ethfunc(B,mf1k(i),alpha1k(i),gJ,gI,Ia,Ehf)/E6 + Ethfunc(B,mf2k(i),alpha2k(i),gJ,gI,Ia,Ehf)/E6;
        chi_array(i) = real(chilfunc(E+Eth0-Eth,ll));
    end
    chi_matrix = diag(chi_array);
    det_value = det(chi_matrix-Kc);
end
