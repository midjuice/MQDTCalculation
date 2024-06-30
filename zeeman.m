clear
clc

element = 'Li7';
channel = 'aa';
ll = 0;
gJ = 2.0023;

[Ia, mf1k, alpha1k, mf2k, alpha2k, gI, sgn, Ehf, E6, beta6, ~] = Get_Params(element, channel);
Bs = linspace(.1,1000,1000);
mf = [1;0;-1;-2;-1;0;1;2];
order = [1;1;1;0;-1;-1;-1;0];

for i = 1:length(Bs) 
    for j = 1:length(mf)
        Eth(i,j)=Ethfunc(Bs(i),mf(j),order(j),gJ,gI,Ia,Ehf)/E6;
    end

end


for k = 1:length(mf)
    plot(Bs,Eth(:,k));
    hold on;
end

ylabel('zeeman splitting energy E_A/E_6');
xlabel('magnet field Bs')