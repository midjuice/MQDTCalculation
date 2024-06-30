function y=Keffc(mf1,alpha1,mf2,alpha2,mf1p,alpha1p,mf2p,alpha2p,Ia1,x1,Ia2,x2,Ks,Kt)
sgnk=[-1 1];

% alpha denotes (M_F)_{1 or 0 or -1}
% beta denotes M_j1 or M_j2
% M_I is independent of K^{c(JI)} matrix, therefore is absent from the
% transform matrix U^{JI}

if(mf1+mf2==mf1p+mf2p && abs(mf1-mf1p)<2 && abs(mf2-mf2p)<2)
    tempKc=0;
    for kk=1:2
        beta1=sgnk(kk);
        for jj=1:2
            beta2=sgnk(jj);
            for ss=1:2
                beta1p=sgnk(ss);
                for tt=1:2
                    beta2p=sgnk(tt);
                    if(2*mf1-beta1==2*mf1p-beta1p && 2*mf2-beta2==2*mf2p-beta2p)
                        % this condition requires that Ufunc1 =
                        % Ufunc3^{\dagger} and Ufunc2 = Ufunc4^{\dagger}
                        temp=Ufunc(alpha1,beta1,mf1,Ia1,x1)*Ufunc(alpha2,beta2,mf2,Ia2,x2)*Ufunc(alpha1p,beta1p,mf1p,Ia1,x1)*Ufunc(alpha2p,beta2p,mf2p,Ia2,x2)*KcJfunc(beta1,beta2,beta1p,beta2p,Ks,Kt);
                        tempKc=tempKc+temp;
                    end
                end
            end            
        end
    end        
    Kc=tempKc;
else
    Kc=0;
end
y=Kc;
end
