clear
clc

% The function below is used to compare the overall inner state energy
% difference between different channels of any nuclide.


element = 'Na23';
channel = 'cc';
gJ = 2.0023;

[Ia, mf1k, alpha1k, mf2k, alpha2k, gI, sgn, Ehf, E6, beta6, a0] = Get_Params(element, channel);
KcST = (0.5*gammacomplex(3/4)+a0./beta6*gammacomplex(5/4))*tan(pi/8)./(-0.5*gammacomplex(3/4)+a0./beta6*gammacomplex(5/4));
Ks = KcST(1);
Kt = KcST(2);

B = linspace(.1, 1000, 1000);
Delta_E_ch1 = zeros(1,1);
Delta_E_ch2 = zeros(1,1);
for i=1:length(B)
    Delta_E_ch1(i) = Ethfunc(B(i), mf2k(1), alpha2k(1), gJ, gI, Ia, Ehf) + Ethfunc(B(i), mf1k(1), alpha1k(1), gJ, gI, Ia, Ehf);
end

for i=2:length(mf1k)
    for j=1:length(B)
        Delta_E_ch2(j) = Ethfunc(B(j), mf2k(i), alpha2k(i), gJ, gI, Ia, Ehf) + Ethfunc(B(j), mf1k(i), alpha1k(i), gJ, gI, Ia, Ehf);
    end
    diff = Delta_E_ch2-Delta_E_ch1;
    plot(B, diff);
    [~, index] = min(diff);
    min_B = B(index);
    xline(min_B, '--', num2str(min_B));
    hold on;
end
