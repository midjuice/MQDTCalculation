clear
clc

element = 'Li7';
channel = 'aa';
ll = 2;
target = 2;
switch target
    case 1
        % B0=1039.24Gs
        left = 1038.8;
        right = 1039.7;
        split_B = [-0.348535345301122,-0.174267672650561,0,0.348535345301122];
    case 2
        % B0=1055.64Gs
        left = 1055.3;
        right = 1056.3;
        split_B = [-0.343780372588848,-0.171890186294424,0,0.343780372588848];
end
points = 100000;
gJ = 2.0023;
[Ia, mf1k, alpha1k, mf2k, alpha2k, gI, sgn, Ehf, E6, mu, C6, a0] = Get_Params(element, channel);
beta6 = (C6*2*mu)^0.25;
KcST = (0.5*gammacomplex(3/4)+a0./beta6*gammacomplex(5/4))*tan(pi/8)./(-0.5*gammacomplex(3/4)+a0./beta6*gammacomplex(5/4));
Ks = KcST(1);
Kt = KcST(2);

Bs = linspace(left, right, points);
output = zeros(1, points);

ppm = ParforProgressbar(points, 'showWorkerProgress', true);
parfor i=1:points
    B = Bs(i);
    output(i) = real((-1)^ll+1/KeffuncID(0,ll,B,mf1k,alpha1k,mf2k,alpha2k,gJ,gI,Ehf,Ia,Ks,Kt,E6,sgn));
    ppm.increment();
end
delete(ppm);

for k=1:length(split_B)
    plot(Bs+split_B(k), output, 'o', markersize=1);
    hold on;
end

xlabel('Magnetic field B/G');
ylabel('Scattering length a/a_0')