clc
clear


element = 'Li7';
channel = 'aa';
ll = 2;
gJ = 2.0023;
[~, ~, ~, ~, ~, ~, ~, ~, E6, mu, C6, ~] = Get_Params(element, channel);
beta6 = (C6*2*mu)^0.25;
channelNum = 5;
B=1039.24;
Nmax = 10;
[Delta,v] = deltaCal(element, channel, ll, B, gJ);
[bjpp_matrix,bjmm_matrix] = bmCal(element, channel, ll, B, Nmax, gJ);
Kc=Kcfunc(element, channel, B, ll, gJ);

points = 10000;
ppm = ParforProgressbar(points, 'showWorkerProgress',true,'parpool', {'local', 6});

r = linspace(.1,2000,points);
parfor i=1:length(r)
    phi(:,i) = phifunc(r(i),element,channel,ll,B,gJ,Nmax,beta6);
    ppm.increment();
end

delete(ppm);

radisplit = 0;
dl = 0.5;
phisq = zeros(channelNum,1);

for i = 1:length(phi)
    index = sum(conj(phi(:,i)).*phi(:,i)*1/r(i));
    radisplit = radisplit + index;
end

split_E = E6*radisplit*anguSplit(ll)*dl;
if B < 1040
    delta_mu = 5.49;
else
    delta_mu = 5.53;
end
split_B = split_E/delta_mu;

disp(split_B);
