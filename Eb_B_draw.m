clear
clc

element = 'Li7';
channel = 'aa';
ll = 2;
gJ = 2.0023;
[Ia, mf1k, alpha1k, mf2k, alpha2k, gI, sgn, Ehf, E6, mu, C6, a0] = Get_Params(element, channel);
beta6 = (C6*2*mu)^0.25;
R_vdw = 1/2*(2*C6*mu)^0.25;
KcST = (0.5*gammacomplex(3/4)+a0./beta6*gammacomplex(5/4))*tan(pi/8)./(-0.5*gammacomplex(3/4)+a0./beta6*gammacomplex(5/4));
Ks = KcST(1);
Kt = KcST(2);

%%%%% This is the main function to calculate the relationship of binding 
%%%%% energy and external applied magnetic field, i.e. E_b-B curve. 
%%%%% "E_range" is the binding energy range to search E_b for fixed B, of
%%%%% which the lower limit can be set as -10e-3 and the upper limit can be
%%%%% set as 10e-3. "B_array" is the magnetic field range near the
%%%%% resonance. "E_value" is used to record the binding energy E_b for a
%%%%% particular B.

points = 100;
left = 1037; %1039.15; %1055.4;
right = 1041; %1039.35; %1056.1;
fit_points = 1000;
E_range = linspace(-0.2, 0.2, 100);
B_array = linspace(left, right, points);
E_value = zeros(1, points);
min_det_value = zeros(1, points);

ppm = ParforProgressbar(points, 'showWorkerProgress', true);
parfor i=1:points
    B = B_array(i);
    [E_value(i), min_det_value(i)] = search_E_value(E_range, ll, sgn, B, mf1k, alpha1k, mf2k, alpha2k, gJ, gI, Ia, Ehf, E6, Ks, Kt)
    ppm.increment();
end
delete(ppm);

plot(B_array, E_value, 'ob', 'MarkerSize', 4, 'MarkerFaceColor', 'b');
hold on;
xline(1039.24,'r--');
yline(0,'r--');


ylabel('Bound state energy E/E_6');
xlabel('Magnetic field B/G');


p = polyfit(B_array, E_value, 1);
Bs = linspace(left, right, fit_points);
E_fit = polyval(p, Bs);

plot(Bs,E_fit,'g','LineWidth',1.5);

[~, ~, ~, ~, ~, ~, ~, ~, E6, ~] = Get_Params("Li7", "aa");

%%%%% CALCULATION OF DELTA_NU, S_RES, R_BG AND ZETA

% \delta\mu in the unit of MHz/Gs
delta_mu_1 = E6*p(1);
% \delta\mu in the unit of Bohr magneton \mu_B, data quoted from
% https://physics.nist.gov/cgi-bin/cuu/Value?mubshhz 
delta_mu_2 = E6*p(1)/1.39962449361;

% a_{bg} and \Delta are fit parameters of the two resonances
abg = 1.529;
delta =  0.0615;

% \bar{a} is the generalized mean scattering length
abar = pi^2/2^(4*ll+1)/(gammacomplex(ll/2+1/4)*gammacomplex(ll+3/2))^2*beta6^(2*ll+1);

% \bar{a}_0 and \bar{E} are in the units of Hartree atomic unit
abar0 = 4*pi*(gammacomplex(1/4)^(-2))*R_vdw;
Ebar = 1/(2*mu*abar0^2);

% r_{bg} is a dimensionless quantity
rbg = abg./abar;

% p(1)*delta is \delta\mu in the units of Hartree atomic unit
sres = rbg.*p(1).*delta./Ebar;

% \zeta is a dimension less quantity
zeta = 1/2.*sres.*abs(rbg);




