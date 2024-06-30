clear
clc

element = 'Li6';
channel = 'ab';
ll = 0;
gJ = 2.0023;

[Ia, mf1k, alpha1k, mf2k, alpha2k, gI, sgn, Ehf, E6, beta6, ~] = Get_Params(element, channel);

% a0 = [44.7634 -2134.8];
% a0 = [42.24 100];

target_B = 834;
N_as = 3;
N_at = 5;
as = linspace(43.6, 44.1, N_as);
at = linspace(-2100, -1100, N_at);
Ks_series = (0.5*gammacomplex(3/4)+as./beta6*gammacomplex(5/4))*tan(pi/8)./(-0.5*gammacomplex(3/4)+as./beta6*gammacomplex(5/4));
Kt_series = (0.5*gammacomplex(3/4)+at./beta6*gammacomplex(5/4))*tan(pi/8)./(-0.5*gammacomplex(3/4)+at./beta6*gammacomplex(5/4));

resonance_loc = zeros(N_at, N_as);

left = .1;
right = 1500;
points = 1000;
Bs = linspace(left, right, points);
output = zeros(1, points);

% ppm = ParforProgressbar(points); 

for m=1:N_as
    for n=1:N_at
        Ks = Ks_series(m);
        Kt = Kt_series(n);
        for i=1:points
            B = Bs(i);
            output(i) = real((-1)^ll+1/KeffuncID(0,ll,B,mf1k,alpha1k,mf2k,alpha2k,gJ,gI,Ehf,Ia,Ks,Kt,E6,sgn));
        end
        [~, Index] = max(output);
        resonance_loc(n, m) = (Bs(Index)+Bs(Index+1))/2;
%         subplot(N_as, N_at, (m-1)*N_as+n);
%         plot(Bs, output, 'ob', markersize=1.5);
%         xline(resonance_loc(n,m), '--', num2str(resonance_loc(n,m)) );
%         title(['a_s = ' num2str(as(m)) ', a_t = ' num2str(at(n))]);
    end
end

% delete(ppm);

surf(as, at, resonance_loc);
hold on;
Z=ones(N_at, N_as)*target_B;
surf(as,at,Z);
xlabel('a_s');
ylabel('a_t');
zlabel('B_{max} location');




