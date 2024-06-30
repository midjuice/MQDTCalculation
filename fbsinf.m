function [finf,ginf] = fbsinf(r,element,channel,ll,B,gJ,Nmax,beta6)
    [Delta,~] = deltaCal(element,channel,ll,B,gJ);
    [Wfp, Wfm, Wgp, Wgm] = Wmatrix(element,channel,ll,B,gJ,Nmax);

    kappa = (-Delta).^(1/2);
    for i = 1:length(Delta)
        finf(i) = (pi* Delta(i))^(-1/2)*(Wfp(i)* exp(-kappa(i)*r/beta6) + Wfm(i)* exp(kappa(i)*r/beta6))/r;
        ginf(i) = (pi* Delta(i))^(-1/2)*(Wgp(i)* exp(-kappa(i)*r/beta6) + Wgm(i)* exp(kappa(i)*r/beta6))/r;
    end
    finf = finf(2:end);
    ginf = ginf(2:end);