function [f,g] = scatterfunc(r,element,channel,ll,B,gJ,Nmax,beta6)

    %%%%    parameter
    [~,v] = deltaCal(element,channel,ll,B,gJ);
    [bjpp_matrix,bjmm_matrix] = bmCal(element,channel,ll,B,Nmax,gJ);

    %%%%    wavefunction
    if r < 10
        [fo,go] = fbarorin(element,channel,ll,B,gJ,Nmax,r);
        f = fo(1);
        g = go(1);
    end
    %%%%    connection 
    if r > 10 && r < 12
        [fi,gi] = fbarmid(r,beta6,v,bjpp_matrix,bjmm_matrix);
        [fo,go] = fbarorin(element,channel,ll,B,gJ,Nmax,r);
        if r > 10.2626
            f = fi(1);
        else
            f = fo(1);
        end

        if r > 11.7374
            g = gi(1);
        else
            g = go(1);
        end
    end
    
    if r > 12
        [fi,gi] = fbarmid(r,beta6,v,bjpp_matrix,bjmm_matrix);
        f =fi(1);
        g =gi(1);
    end
            
end


