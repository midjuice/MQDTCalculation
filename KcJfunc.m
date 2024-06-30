function y=KcJfunc(beta1,beta2,beta1p,beta2p,Ks,Kt)
temp1=cgfuncST(beta1,beta2,0,0)*cgfuncST(beta1p,beta2p,0,0);
temp2=0;
for mj=-1:1
    temp2=temp2+cgfuncST(beta1,beta2,1,mj)*cgfuncST(beta1p,beta2p,1,mj);
end
y=Ks*temp1+Kt*temp2;
end