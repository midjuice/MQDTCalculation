clc
clear

% fail attempt to draw binding energy spectrum
% We used a wrong equation E_b=-\hbar^2/(2*m*a^2), which only applies to 
% the case of s wave resonance

element = 'Li7';
channel = 'aa';
ll = 2;
gJ = 2.0023;

left2 = 1069.08;
right2 = 1069.12;
points = 1000;
left1 = 1052.97;
right1 = 1052.99;


[Ia, mf1k, alpha1k, mf2k, alpha2k, gI, sgn, Ehf, E6, mu, C6, a0] = Get_Params(element, channel);
beta6 = (C6*2*mu)^0.25;
R_vdw = 1/2*(2*C6*mu)^0.25;

KcST = (0.5*gammacomplex(3/4)+a0./beta6*gammacomplex(5/4))*tan(pi/8)./(-0.5*gammacomplex(3/4)+a0./beta6*gammacomplex(5/4));
Ks = KcST(1);
Kt = KcST(2);
Bs1 = linspace(left1, right1, points);
Bs2 = linspace(left2, right2, points);
a1 = zeros(1, points);
a2 = zeros(1, points);


ppm = ParforProgressbar(points, 'showWorkerProgress',true,'parpool', {'local', 6});
parfor i=1:points
        B = Bs1(i);
        a1(i) = real((-1)^ll+1/KeffuncID(0,ll,B,mf1k,alpha1k,mf2k,alpha2k,gJ,gI,Ehf,Ia,Ks,Kt,E6,sgn));
        ppm.increment();
end

parfor i=1:points
        B = Bs2(i);
        a2(i) = real((-1)^ll+1/KeffuncID(0,ll,B,mf1k,alpha1k,mf2k,alpha2k,gJ,gI,Ehf,Ia,Ks,Kt,E6,sgn));

        ppm.increment();
end
[~, Index1] = max(a1);
[~, Index2] = max(a2);
Eb1 = -1./(2*mu.*a1.^2);
Eb2 = -1./(2*mu.*a2.^2);


delete(ppm);

p1 = polyfit(Bs1(Index1-40:Index1-10),Eb1(Index1-40:Index1-10),1);
x1 = linspace(Bs1(Index1-40),Bs1(Index1-10),100);
y1 = polyval(p1,x1); 

p2 = polyfit(Bs2(Index2-100:Index2-50),Eb2(Index2-100:Index2-50),1);
x2 = linspace(Bs2(Index2-100),Bs2(Index2-20),100);
y2 = polyval(p2,x2); 




%%%%% PLOT OF BINDING ENERGY AND FITTING
tiledlayout(2,2);

nexttile;
plot(x1,y1,'g',markersize = 8);
hold on;
plot(Bs1(Index1-40:Index1),Eb1(Index1-40:Index1),'ob',markerSize = 1.5);
title('binding energy');
ylabel('binding energy');


nexttile;
plot(x2,y2,'g',markersize = 8);
hold on;
plot(Bs2(Index2-100:Index2),Eb2(Index2-100:Index2),'ob',markerSize = 1.5);

title('binding energy');
ylabel('binding energy');


nexttile;
plot(Bs1,a1,'or',markerSize = 1.5);

nexttile;
plot(Bs2,a2,'or',markerSize = 1.5);
title('scattering length');

xlabel('B');
ylabel('scattering length');

%%%%% FITTING  RESULT OF ABG AND DELTA
abg = [1.536,1.53];
delta = [0.0006223, 0.05871];

%%%%% CALCULATION OF DELTA_NU, S_RES, R_BG AND ZETA
deltanu = [p1(1),p2(1)];

% corresponding energy scale
abar = 4*pi*(gamma(1/4)^(-2))*R_vdw;
Ebar = 1/(2*mu*abar^2);

rbg = abg./abar;
sres = rbg.*deltanu.*delta./Ebar;
zeta = 1/2.*sres.*abs(rbg);




