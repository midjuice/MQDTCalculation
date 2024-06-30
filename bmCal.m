function [bjpp_matrix,bjmm_matrix] = bmCal(element,channel,ll,B,Nmax,gJ)
    
    [Delta,v] = deltaCal(element,channel,ll,B,gJ);

    Nb = zeros(1,1);
    bjpp_matrix = zeros(5,18);
    bjmm_matrix = zeros(5,18);
    for i =1:length(v)
        [bjpp,bjmm,Nb(i),~,~]=bjfunc(v(i),Delta(i),ll,Nmax);
        bjpp_matrix(i,1:length(bjpp)) = bjpp;
        bjmm_matrix(i,1:length(bjmm)) = bjmm;
    end
end