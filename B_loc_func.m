function y = B_loc_func(element, channel_series, ll_series, left_series, right_series, points_series, gJ, opt_params)
    num_resonance = length(channel_series);
    y = zeros(1, num_resonance);
    a0 = zeros(1, 2);
    a0(1) = opt_params(1);
    a0(2) = opt_params(2);
    C6 = opt_params(3);
    parfor s=1:num_resonance
        channel = channel_series(s);
        ll = ll_series(s);
        left = left_series(s);
        right = right_series(s);
        points = points_series(s);

        [Ia, mf1k, alpha1k, mf2k, alpha2k, gI, sgn, Ehf, E6, mu, ~] = Get_Params(element, channel);
        beta6 = (C6*2*mu)^0.25;
        KcST = (0.5*gammacomplex(3/4)+a0./beta6*gammacomplex(5/4))*tan(pi/8)./(-0.5*gammacomplex(3/4)+a0./beta6*gammacomplex(5/4));
        Ks = KcST(1);
        Kt = KcST(2);
        Bs = linspace(left, right, points);
        output = zeros(1, points);
        
        for i=1:points
            B = Bs(i);
            output(i) = real((-1)^ll+1/KeffuncID(0,ll,B,mf1k,alpha1k,mf2k,alpha2k,gJ,gI,Ehf,Ia,Ks,Kt,E6,sgn));
        end
        [~, Index] = max(output);
        y(s) = (Bs(Index)+Bs(Index+1))/2;
    end
end