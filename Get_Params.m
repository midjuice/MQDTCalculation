function [Ia, mf1k, alpha1k, mf2k, alpha2k, gI, sgn, Ehf, E6, mu, C6, a0] = Get_Params(element, channel)

switch element
    case {'Na23', 'Li7', 'K39', 'Rb87', 'K41'}
        Ia = 3/2;
        switch channel
            case 'aa'
                mf1k=   [ 1   0   0   1   1];
                alpha1k=[ 1   1  -1   1  -1];
                mf2k=   [ 1   2   2   1   1];
                alpha2k=[ 1   0   0  -1  -1];
            case 'cc'
                mf1k=   [-1   0   0  -1  -1];
                alpha1k=[ 1   1  -1   1  -1];
                mf2k=   [-1  -2  -2  -1  -1];
                alpha2k=[ 1   0   0  -1  -1];
            case 'ae'
                mf1k   =[  1     1    0    2 ];
                alpha1k=[  1    -1    -1   0 ];
                mf2k=   [  -1   -1     0  -2 ];
                alpha2k=[  -1   -1    -1   0 ];
            case 'ab'
                mf1k   =[  1     1    1 ];
                alpha1k=[  1    1     -1];
                mf2k=   [  0    0     0];
                alpha2k=[  1    -1    -1];
            case 'ac'
                mf1k=   [ 1   0   0  -1  -1];
                alpha1k=[ 1   1  -1   1  -1];
                mf2k=   [-1  -2  -2  -1  -1];
                alpha2k=[ 1   0   0  -1  -1];
            case 'bc'
                mf1k=   [ 0   0   0   0  ];
                alpha1k=[ 1   1  -1  -1  ];
                mf2k=   [-1  -1  -1  -1  ];
                alpha2k=[ 1  -1   1  -1  ];
        end

    case 'Li6'
        Ia = 1;
        switch channel
            case 'aa'
                mf1k   =[1/2  1/2  1/2  3/2   -1/2 ];
                alpha1k=[ 1,   1,   -1,  0,   -1   ];
                mf2k=   [1/2  1/2  1/2  -1/2   3/2 ];
                alpha2k=[ 1,   -1,  -1,  1,    0   ];
            case 'ab'
                mf1k    = [1/2   1/2  1/2  1/2  3/2];
                alpha1k = [ 1,   1,   -1,  -1,   0 ];
                mf2k    = [-1/2 -1/2 -1/2 -1/2 -3/2];
                alpha2k = [ 1,  -1,    1,  -1,   0 ];
            case 'ac'
                mf1k    = [ 1/2   1/2   -1/2  -1/2];
                alpha1k = [ 1     -1     1      1];
                mf2k    = [-3/2  -3/2   -1/2   -1/2];
                alpha2k = [ 0      0    -1     1];
            case 'bb'
                mf1k    = [-1/2   1/2    1/2  -1/2  -1/2];
                alpha1k = [ 1     -1     1    -1    -1  ];
                mf2k    = [-1/2  -3/2   -3/2  -1/2  -1/2];
                alpha2k = [ 1      0     0    -1     1  ];
            case 'bc'
                mf1k    = [-1/2  -1/2];
                alpha1k = [ 1     -1 ];
                mf2k    = [-3/2  -3/2];
                alpha2k = [ 0      0 ];
        end
    case 'K40'
        Ia = 4;
        switch channel
            case 'ab'
                mf1k   =[  -9/2    -9/2      ];
                alpha1k=[   0        0       ];
                mf2k=   [  -7/2    -7/2      ];
                alpha2k=[  -1        1       ];
            case 'ac'
                mf1k   =[  -5/2    -9/2    -7/2    -7/2];
                alpha1k=[  -1        0       1        1];
                mf2k=   [  -9/2    -5/2     -7/2   -7/2];
                alpha2k=[   0        1       -1       1];
            case 'bb'
                mf1k   =[  -7/2    -7/2   -7/2    -5/2   -5/2];
                alpha1k=[  -1        1     -1      1     -1];
                mf2k=   [  -7/2    -7/2   -7/2    -9/2   -9/2];
                alpha2k=[  -1        1      1       0      0];
        end
    case 'Rb85'
        Ia = 5/2;
        switch channel
            case 'ee'
                mf1k   =[ -2    -2   -2     -1   -1];
                alpha1k=[  1    -1   -1     1   -1];
                mf2k=   [ -2    -2   -2    -3   -3];
                alpha2k=[  1    -1   1      0    0];
        end
    case 'Cs133'
        Ia = 7/2;
        switch channel
            case 'aa'
                mf1k   =[  -3  -3    -3    -4   -4   ];
                alpha1k=[   1   -1    -1    0    0   ];
                mf2k=   [  -3  -3    -3    -2   -2   ];
                alpha2k=[   1   1    -1    -1    1   ];
        end
end

switch element
    case 'Li6'
        gI=-4.477e-4;
        sgn=-1;
        Ehf=228.205;
        m=6.01512;
        C6 = 1039.39;%1404.37614014271; %1393.39;
        a0 = [44.7634 -2134.8];%[44.6030599346324	-2041.06198087401]; 
    case 'Li7'
        gI=-5.429e-4;
        sgn=1;
        Ehf=803.504;
        m=7.01601;
        C6 = 1393.39;% 1446.76482808662; %1393.39;
        a0 = [34.331 -26.92];%[33.7368999427788	-27.8960181429473]; %[34.331 -26.92];
    case 'Na23'
        gI=-8.0466e-04;
        sgn=1;
        Ehf=1772;
        m=22.9898;
        C6 = 1556;
        a0 = [15.5 62.7];
    case 'K39'
        gI=-1.419e-04;
        sgn=1;
        Ehf=230.8599;
        m=38.963707;
        C6 = 3825.6;
        a0 = [177 -15.5];
    case 'K40'
        gI = 0.000176490;
        sgn = -1;
        Ehf = -285.7308;
        m = 39.96399;
        C6 = 3897;
        a0 = [-885 204];
    case 'K41'
        gI = -7.7913e-05;
        sgn = 1;
        Ehf = 127.0069352;
        m = 40.961825;
        C6 = 3897;
        a0 = [82.35 60.305];
    case 'Rb85'
        gI=-2.9370e-04;
        sgn=1;
        Ehf=1011.91;
        m = 84.9118;
        C6 = 4698;
        a0 = [2400 -369];
    case 'Rb87'
        gI=-9.9535e-04;
        sgn=1;
        Ehf=3417.34;
        m = 86.909184;
        C6 = 4698;
        a0 = [90 106];
    case 'Cs133'
        gI=-0.000398;
        sgn=1;
        Ehf=2298;
        m=132.90545193;
        C6 = 6860;
        a0 = [280.37 2440];
end

mu=m/2*1822.89;
R_vdw = 1/2*(2*C6*mu)^0.25;
E_vdw = 1/(2*mu)/R_vdw^2*6579.684*1e6;
E6 = E_vdw/4;

end

%%%%% "C6", "mu", "R_vdw" are in Hartree atomic units. "E_vdw", "E6" are in
%%%%% the units of mega-Hertz.