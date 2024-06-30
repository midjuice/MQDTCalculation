clear
clc

element = 'Li7';
channel = 'aa';
ll = 0;
left = 730;
right = 750;
points = 1000;

gJ = 2.0023;
[Ia, mf1k, alpha1k, mf2k, alpha2k, gI, sgn, Ehf, E6, mu, C6, a0] = Get_Params(element, channel);
beta6 = (C6*2*mu)^0.25;
KcST = (0.5*gammacomplex(3/4)+a0./beta6*gammacomplex(5/4))*tan(pi/8)./(-0.5*gammacomplex(3/4)+a0./beta6*gammacomplex(5/4));
Ks = KcST(1);
Kt = KcST(2);
R_vdw = 1/2*(2*C6*mu)^0.25;
beta6 = (C6*2*mu)^0.25;
E = zeros(1,points);

Bs = linspace(left, right, points);
output = zeros(1, points);

ppm = ParforProgressbar(points, 'showWorkerProgress', true);
parfor i=1:points
    B = Bs(i);
    output(i) = ebfunc(ll,B,mf1k,alpha1k,mf2k,alpha2k,gJ,gI,Ehf,Ia,Ks,Kt,E6,sgn)
    ppm.increment();
end
delete(ppm);
plot(Bs, output, 'ob', markersize=1.5);




