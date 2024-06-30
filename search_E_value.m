function [E_value, min_det_value] = search_E_value(E_range, ll, sgn, B, mf1k, alpha1k, mf2k, alpha2k, gJ, gI, Ia, Ehf, E6, Ks, Kt)
    det_value = zeros(1, length(E_range));
    for i=1:length(E_range)
         det_value(i) = bound_spectrum_eqn(E_range(i), ll, sgn, B, mf1k, alpha1k, mf2k, alpha2k, gJ, gI, Ia, Ehf, E6, Ks, Kt);
    end
    [~, min_det_index] = min(abs(det_value));
    min_det_value = det_value(min_det_index);
    E_value = E_range(min_det_index);
end