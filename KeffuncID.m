function y=KeffuncID(Es,ll,B,mf1k,alpha1k,mf2k,alpha2k,gJ,gI,Ehf,Ia,Ks,Kt,E6,sgn,~,~)
Kc=zeros(1,1);
x=(gJ-gI)*B*1.399624504/Ehf;

Nc=max(size(mf1k));
for kk=1:Nc
    for jj=1:Nc
        temp1=Keffc(mf1k(kk),alpha1k(kk),mf2k(kk),alpha2k(kk),mf1k(jj),alpha1k(jj),mf2k(jj),alpha2k(jj),Ia,x,Ia,x,Ks,Kt);
        temp2=Keffc(mf2k(kk),alpha2k(kk),mf1k(kk),alpha1k(kk),mf1k(jj),alpha1k(jj),mf2k(jj),alpha2k(jj),Ia,x,Ia,x,Ks,Kt);
        temp3=Keffc(mf1k(kk),alpha1k(kk),mf2k(kk),alpha2k(kk),mf2k(jj),alpha2k(jj),mf1k(jj),alpha1k(jj),Ia,x,Ia,x,Ks,Kt);
        temp4=Keffc(mf2k(kk),alpha2k(kk),mf1k(kk),alpha1k(kk),mf2k(jj),alpha2k(jj),mf1k(jj),alpha1k(jj),Ia,x,Ia,x,Ks,Kt);
        Kc(kk,jj)=(temp1+temp4+sgn*(-1)^ll*(temp3+temp2))/2;
        if(mf1k(kk)==mf2k(kk) && alpha1k(kk)==alpha2k(kk))
            Kc(kk,jj)=Kc(kk,jj)/sqrt(2);
        end
        if(mf1k(jj)==mf2k(jj) && alpha1k(jj)==alpha2k(jj))
            Kc(kk,jj)=Kc(kk,jj)/sqrt(2);
        end
    end
end
%Kc

Eth10=Ethfunc(B,mf1k(1),alpha1k(1),gJ,gI,Ia,Ehf)/E6;
Eth20=Ethfunc(B,mf2k(1),alpha2k(1),gJ,gI,Ia,Ehf)/E6;
Eth0=Eth10+Eth20;


if(Nc>1)
chik=zeros(1,1);
K1c=zeros(1,1);
for kk=2:Nc
    Eth1k=Ethfunc(B,mf1k(kk),alpha1k(kk),gJ,gI,Ia,Ehf)/E6;
    Eth2k=Ethfunc(B,mf2k(kk),alpha2k(kk),gJ,gI,Ia,Ehf)/E6;
    Ethk=Eth1k+Eth2k;%Ethk*E6
    chik(kk-1,kk-1)=chilfunc(Es+Eth0-Ethk,ll);
    K1c(1,kk-1)=Kc(1,kk);
end

Kcc=zeros(1,1);
for kk=2:Nc
    for jj=2:Nc
        Kcc(kk-1,jj-1)=Kc(kk,jj);
    end
end

Keffc0=Kc(1,1)+(K1c/(chik-Kcc))*(K1c');
else
    Keffc0=Kc(1,1);
end

nu0=(2*ll+1)/4;
y=(Keffc0-tan(pi*nu0/2))/(1+tan(pi*nu0/2)*Keffc0);
end