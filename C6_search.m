clear
clc

element = 'Li7';
channel = 'aa';
ll = 0;
gJ = 2.0023;

[Ia, mf1k, alpha1k, mf2k, alpha2k, gI, sgn, Ehf, E6, ~] = Get_Params(element, channel);


a0 = [34.331 -26.92];
as = a0(1);
at = a0(2);
beta6 = linspace(50, 70, 100);
% a0 = [42.24 100];

target_B = 736;
target_line = ones(1, length(beta6))*target_B;

Ks_series = (0.5*gammacomplex(3/4)+as./beta6.*gammacomplex(5/4)).*tan(pi/8)./(-0.5*gammacomplex(3/4)+as./beta6.*gammacomplex(5/4));
Kt_series = (0.5*gammacomplex(3/4)+at./beta6.*gammacomplex(5/4)).*tan(pi/8)./(-0.5*gammacomplex(3/4)+at./beta6.*gammacomplex(5/4));

resonance_loc = zeros(1, length(beta6));

left = 600;
right = 900;
points = 1000;
Bs = linspace(left, right, points);
output = zeros(1, points);

% ppm = ParforProgressbar(points); 

for n=1:length(beta6)
    Ks = Ks_series(n);
    Kt = Kt_series(n);
    parfor i=1:points
        B = Bs(i);
        output(i) = real((-1)^ll+1/KeffuncID(0,ll,B,mf1k,alpha1k,mf2k,alpha2k,gJ,gI,Ehf,Ia,Ks,Kt,E6,sgn));
    end
    [~, Index] = max(output);
    resonance_loc(1, n) = Bs(Index);
%         subplot(N_as, N_at, (m-1)*N_as+n);
%         plot(Bs, output, 'ob', markersize=1.5);
%         xline(resonance_loc(n,m), '--', num2str(resonance_loc(n,m)) );
%         title(['a_s = ' num2str(as(m)) ', a_t = ' num2str(at(n))]);
end


% delete(ppm);

plot(beta6, resonance_loc, 'o-');
hold on;
plot(beta6, target_line);




