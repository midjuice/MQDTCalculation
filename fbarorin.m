function [f,g] = fbarorin(element,channel,ll,B,gJ,Nmax,r)


    [alphal,betal] = abCal(element,channel,ll,B,gJ,Nmax);
    vo=(2*ll+1)/4;
    [~, ~, ~, ~, ~, ~, ~, ~, ~, mu, C6, ~] = Get_Params(element, channel);
    beta6 = (C6*2*mu)^0.25;

    f = sqrt(4/pi) * (r/beta6) * r^(-1/2) * (alphal* cos( 0.5*(r/beta6)^(-2)-vo*pi/2-pi/4 )  + betal * sin( 0.5*(r/beta6)^(-2)-vo*pi/2-pi/4 ) );
    g = sqrt(4/pi) * (r/beta6) * r^(-1/2) * (-betal* cos( 0.5*(r/beta6)^(-2)-vo*pi/2-pi/4 )  + alphal * sin( 0.5*(r/beta6)^(-2)-vo*pi/2-pi/4 ) );
end