function res = anguSplit(ll)
    switch ll
        case 0
            res = 0;
        case 1
            res = [-4/5,2/5];
        case 2 
            res = [-4/7,-2/7,4/7];
    end
end


%%%%% result from mathematica, which code is Integrate[SphericalHarmonicY[1, 1, \[Theta], \[Phi]]*SphericalHarmonicY[1, 1, \[Theta], -\[Phi]]*
%%%%% Sin[\[Theta]]*(1 - 3*Cos[\[Theta]]^2), {\[Theta], 0, Pi}, {\[Phi], 0, 2*Pi}]