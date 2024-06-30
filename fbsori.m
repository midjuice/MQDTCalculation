function [fori,gori] = fbsori(r,element,channel,gJ,ll,B,beta6)
    [~,v] = deltaCal(element,channel,ll,B,gJ);
    %%%% parameter
    y = 1/2*(r/beta6)^(3/2);
    vo=(2*ll+1)/4;
    %%%% calculation
    for i = 1:length(v)
        fc0 = (2/pi)^(1/2)*(r/beta6)^(3/2)*cos(y-pi*vo/2-pi/4);
        gc0 = -(2/pi)^(1/2)*(r/beta6)^(3/2)*sin(y-pi*vo/2-pi/4);
        fori(i) = (cos(pi*v(i)/2)*fc0 + sin(pi*v(i)/2)*gc0)/r;
        gori(i) = (-sin(pi*v(i)/2)*fc0 + cos(pi*v(i)/2)*gc0)/r;
    end

    fori = fori(2:end);
    gori = gori(2:end);
end