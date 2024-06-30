clear
clc

element = 'Li7';
channel = 'aa';
ll = 0;
left = 500;
right = 1000;
points = 10000;

fit_judge = 0;
ope_points = 100;
abg_start = 100;
Delta_start = -100;

sub_func_on = 0;
sub_left = 195;
sub_right = 200;
sub_points = 5000;

postion_on = 1;
deltamu_on = 0;
llf = 976;
rrf = 977;

gJ = 2.0023;
[Ia, mf1k, alpha1k, mf2k, alpha2k, gI, sgn, Ehf, E6, mu, C6, a0] = Get_Params(element, channel);
beta6 = (C6*2*mu)^0.25;
KcST = (0.5*gammacomplex(3/4)+a0./beta6*gammacomplex(5/4))*tan(pi/8)./(-0.5*gammacomplex(3/4)+a0./beta6*gammacomplex(5/4));
Ks = KcST(1);
Kt = KcST(2);
R_vdw = 1/2*(2*C6*mu)^0.25;
beta6 = (C6*2*mu)^0.25;

Bs = linspace(left, right, points);
output = zeros(1, points);

ppm = ParforProgressbar(points, 'showWorkerProgress', true);
parfor i=1:points
    B = Bs(i);
    output(i) = real((-1)^ll+1/KeffuncID(0,ll,B,mf1k,alpha1k,mf2k,alpha2k,gJ,gI,Ehf,Ia,Ks,Kt,E6,sgn));
    ppm.increment();
end
delete(ppm);
plot(Bs, output, 'ob', markersize=1.5);
hold on;

%%%%% FUNCTION 1 TO DRAW X-LINES
% peakLoc = islocalmax(output);
% po = Bs(peakLoc);
% t = size(po);
% for i=1:t(2)
%     xline(po(i),'--',po(i));
% end

%%%%% FUNCTION 2 TO DRAW X-LINES
% avg = mean(output);
% temp = output-avg;
% for i=2:points
%     if ( abs(temp(i)-temp(i-1))>1 && temp(i)*temp(i-1)<0 )
%         B_loc = (Bs(i)+Bs(i-1))/2;
%         xline(B_loc,'--',num2str(B_loc));
%     end
% end
% hold on;

%%%%% FUNCTION 3 TO DRAW X-LINES
if (postion_on == 1)
    temp = mean(output);
for i=(ope_points/2+1):(points-ope_points/2)
    if ( abs(output(i)-output(i+1))>5*abs(output(i+10)-output(i-10)) && (output(i)-temp)*(output(i+1)-temp)<0 )
        B_loc = (Bs(i)+Bs(i+1))/2;
        xline(B_loc,'--',num2str(B_loc));
        if fit_judge == 1
            Bs_selected = Bs(1,i-ope_points/2:i+ope_points/2);
            output_selected = output(1,i-ope_points/2:i+ope_points/2);
            f = fittype("abg*(1+delta/(t0-t))", independent='t', coefficients=["abg" "delta" "t0"]);
            fit_fun=fit(Bs_selected', output_selected', f, 'StartPoint', [abg_start, Delta_start, B_loc]);
            disp(fit_fun);
            fitted_output = fit_fun(Bs_selected);
            plot(Bs_selected, fitted_output, 'g');
        end
    end
end
end
hold on;




if ( sub_func_on==1 )
    Bss = linspace(sub_left, sub_right, sub_points);
    outputt = zeros(1,1);
    t = size(Bss);
    for i=1:t(2)
        B = Bss(i);
        outputt(i) = real((-1)^ll+1/KeffuncID(0,ll,B,mf1k,alpha1k,mf2k,alpha2k,gJ,gI,Ehf,Ia,Ks,Kt,E6,sgn));
        if mod(i,sub_points/10)==0
            disp(['Sub-function has finished ', num2str(i/sub_points*100), '%']);
        end
    end
    plot(Bss, outputt, 'ob', markersize=1.5);
    hold on;
    [~, Index] = max(outputt);
    B_loc = (Bss(Index)+Bss(Index+1))/2;
    xline(B_loc,'--',num2str(B_loc));
%     peakLoc = islocalmax(outputt);
%     po = Bss(peakLoc);
%     t = size(po);
%     for i=1:t(2)
%         xline(po(i),'--',po(i));
%     end
end


if (deltamu_on == 1)
    Bde = linspace(llf,rrf, 100);
    for i = 1:length(Bde)
        B = Bde(i);
        output(i) = real((-1)^ll+1/KeffuncID(0,ll,B,mf1k,alpha1k,mf2k,alpha2k,gJ,gI,Ehf,Ia,Ks,Kt,E6,sgn));
        Eb(i) = 1/2/mu/output(i)^2;   
    end

    p = polyfit(Bde, Eb, 1);
    E_fit = polyval(p, Bde);
    % \delta\mu in the unit of MHz/Gs
    delta_mu_1 = E6*p(1);
    % \delta\mu in the unit of Bohr magneton \mu_B, data quoted from
    % https://physics.nist.gov/cgi-bin/cuu/Value?mubshhz 
    delta_mu_2 = E6*p(1)/1.39962449361;

    % a_{bg} and \Delta are fit parameters of the two resonances
    abg = 1.1;
    delta = 121.1;

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
end
% 
% fig = gcf;
% fig.Position(3) = 8;  % 宽度
% fig.Position(4) = 4;  % 高度
xlabel('Magnetic field B/G');
ylabel('Scattering Length a/a_0');
% text(1038.9,10,'A.');

