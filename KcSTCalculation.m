% a0 = [52.50 63.70];
% m1 = 6.01512;
% m2 = 39.9640;
% C6 = 6.01512; % data of Li6-K40, quoted from 10.1103/PhysRevA.89.062718

% a0 = [19.69 64.57];
% m1 = 22.9898;
% m2 = 22.9898;
% C6 = 1556; % data of Na23-Na23, quoted from 10.1103/PhysRevA.72.042719

% a0 = [44.41531297, -1970.318493]; %[44.7634 -2134.8]; %[45.154 -2133];
% m1 = 6.01512;
% m2 = 6.01512;
% C6 = 1393.39; % data of Li6-Li6, quoted from 10.1103/PhysRevA.89.052715

a0 = [34.331 -26.92];
m1 = 7.01601;
m2 = 7.01601;
C6 = 1393.39; % data of Li7-Li7, quoted from 10.1103/PhysRevA.89.052715

% a0 = [278 81.1];
% m1 = 38.963707;
% m2 = 38.963707;
% C6 = 3897; % data of K39-K39, quoted from https://journals.aps.org/pra/pdf/10.1103/PhysRevA.57.R4118

% a0 = [158 1.7];
% m1 = 39.96399;
% m2 = 39.96399;
% C6 = 3897; % data of K40-K40, quoted from https://journals.aps.org/pra/pdf/10.1103/PhysRevA.57.R4118

% a0 = [82.35 59.35]; % [121 286]
% m1 = 40.961825;
% m2 = 40.961825;
% C6 = 3897; % data of K41-K41, quoted from https://journals.aps.org/pra/pdf/10.1103/PhysRevA.57.R4118

% a0 = [280.37 2440];
% m1 = 132.90545193;
% m2 = 132.90545193;
% C6 = 6860; % data of Cs133-Cs133, quoted from 10.1103/PhysRevA.70.032701

mu = m1*m2/(m1+m2)*1822.89;

beta6 = (C6*2*mu)^0.25;
Kcs = ( 0.5*gammacomplex(3/4)+a0./beta6*gammacomplex(5/4) )*tan(pi/8) ./ ( -0.5*gammacomplex(3/4)+a0./beta6*gammacomplex(5/4) );
disp(Kcs)