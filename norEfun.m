function [psi,sumPer] = norEfun(Es,element,channel,B,ll,gJ)
rsize = 1000;
rs = linspace(0.01,100,rsize);
psi = zeros(5,rsize);

for i = 1:length(rs)
    psi(:,i) = ZeroEfun(Es,element,channel,B,ll,gJ,rs(i));
end

psiSq = psi * psi';
sumPsi = sum(psiSq,2);

sumTot = sum(sumPsi);
sumPer = sumPsi/sumTot;



end


