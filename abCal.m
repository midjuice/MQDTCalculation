function [alphal,betal] = abCal(element,channel,ll,B,gJ,Nmax)
    [Delta,v] = deltaCal(element,channel,ll,B,gJ);

    
    alphal = zeros(1,length(v));
    betal = zeros(1,length(v));
    v0=(2*ll+1)/4;
    for i=1:length(v)
        [Xl,Yl,~,~,~] = XYMfunc(v(i),Delta(i),ll,Nmax);
        alphal(i) = cos(pi*(v(i)-v0))*Xl - sin(pi*(v(i)-v0))*Yl;
        betal(i) = sin(pi*(v(i)-v0))*Xl + cos(pi*(v(i)-v0))*Yl;
    end
