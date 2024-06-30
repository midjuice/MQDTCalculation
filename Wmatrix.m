function [Wfp, Wfm, Wgp, Wgm] = Wmatrix(element,channel,ll,B,gJ,Nmax)
    [Delta,v] = deltaCal(element,channel,ll,B,gJ);
    Wfp = zeros(1,length(v));
    Wfm = zeros(1,length(v));
    Wgp = zeros(1,length(v));
    Wgm = zeros(1,length(v));
    for i = 1:length(v)
        [Xl,Yl,Ml,~,Gv] = XYMfunc(v(i),Delta(i),ll,Nmax);
        Wfp(i) = -((Xl^2+Yl^2)*sin(pi*v(i)))^(-1)*Gv*cos(pi*v(i))*((1-Ml)*sin(pi*v(i)/2)*Xl+(1+Ml)*cos(pi*v(i)/2)*Yl);
        Wfm(i) = ((Xl^2+Yl^2)*sin(pi*v(i))*2)^(-1)*Gv/sin(pi*v(i))*((1+Ml)*sin(pi*v(i)/2)*Xl+(1-Ml)*cos(pi*v(i)/2)*Yl);
        Wgp(i) = -((Xl^2+Yl^2)*sin(pi*v(i)))^(-1)*Gv*cos(pi*v(i))*((1+Ml)*cos(pi*v(i)/2)*Xl-(1-Ml)*sin(pi*v(i)/2)*Yl);
        Wgm(i) = ((Xl^2+Yl^2)*sin(pi*v(i))*2)^(-1)*Gv/sin(pi*v(i))*((1-Ml)*cos(pi*v(i)/2)*Xl-(1+Ml)*sin(pi*v(i)/2)*Yl);
    end
