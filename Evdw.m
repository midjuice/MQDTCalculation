m = 40.961825;
C6 = 3897;

mu = m/2*1822.89;
R_vdw = 1/2*(2*C6*mu)^0.25;
E_vdw = 1/(2*mu)/R_vdw^2*6579.684*1e6;
disp(R_vdw);
disp(E_vdw/4);