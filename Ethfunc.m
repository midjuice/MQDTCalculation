function y=Ethfunc(B,mf,alpha,gJ,gI,Ia,Ehf)
x=(gJ-gI)*B*1.399624504/Ehf;
if(abs(2*mf)==2*Ia+1)
    y=(Ia/(2*Ia+1)+sign(mf)*0.5*((gJ+2*Ia*gI)/(gJ-gI))*x)*Ehf;
else
    y=0.5*Ehf*(-1/(2*Ia+1)+2*mf*gI*x/(gJ-gI)-alpha*sqrt(1+4*mf*x/(2*Ia+1)+x^2));
end

end