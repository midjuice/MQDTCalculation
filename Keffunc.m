function y=Keffunc(Es,ll,B,mf1k,alpha1k,mf2k,alpha2k,gJ,gI1,gI2,Ehf1,Ehf2,Ia1,Ia2,Ks,Kt,E6,~,~)
Kc=[0 0;0 0];
x1=(gJ-gI1)*B*1.399624504/Ehf1;
x2=(gJ-gI2)*B*1.399624504/Ehf2;
Nc=max(size(mf1k));
for kk=1:Nc
    for jj=1:Nc
        Kc(kk,jj)=Keffc(mf1k(kk),alpha1k(kk),mf2k(kk),alpha2k(kk),mf1k(jj),alpha1k(jj),mf2k(jj),alpha2k(jj),Ia1,x1,Ia2,x2,Ks,Kt);
    end
end

Eth10=Ethfunc(B,mf1k(1),alpha1k(1),gJ,gI1,Ia1,Ehf1)/E6;
Eth20=Ethfunc(B,mf2k(1),alpha2k(1),gJ,gI2,Ia2,Ehf2)/E6;
Eth0=Eth10+Eth20;


if(Nc>1)
    chik=zeros(1,1);K1c=zeros(1,1);Kcc=zeros(1,1);
    for kk=2:Nc
        Eth1k=Ethfunc(B,mf1k(kk),alpha1k(kk),gJ,gI1,Ia1,Ehf1)/E6;
        Eth2k=Ethfunc(B,mf2k(kk),alpha2k(kk),gJ,gI2,Ia2,Ehf2)/E6;
        Ethk=Eth1k+Eth2k;
        chik(kk-1,kk-1)=chilfunc(Es+Eth0-Ethk,ll);
        K1c(1,kk-1)=Kc(1,kk);
    end
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