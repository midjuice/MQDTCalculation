function dist = dist_func(opt_params, weight)
    element = "Li7";
    channel_series = ["aa"];
    ll_series = [0];
    left_series = [600];
    right_series = [1000];
    points_series = [1000];
    B_loc_target = [737.69];
    gJ = 2.0023;
    B_loc = B_loc_func(element, channel_series, ll_series, left_series, right_series, points_series, gJ, opt_params);
    dist = sum((B_loc_target-B_loc).^2.*weight);
end

% Li7
% B_loc_target = [737.69];
% channel_series = ["aa"];
% ll_series = [0];
% B_loc = [738.223];
% Delta = [-159.8];
% a_{bg} = [-0.6757];


% Li6
% B_loc_target = [832.18, 543.286, 225.33, 689.68, 159.097, 185.091, 214.825, 809.76];
% channel_series = ["ab", "ab", "ac", "ac", "aa", "ab", "bb", "bc"];
% ll_series = [0, 0, 1, 0, 1, 1, 1, 0];
% B_loc = [828.45, 549.70, 206.75, 691.05, 156.503, 182.509, 214.825, 808.16];
% Delta = [-254.7,  0.099, -22.22, -112.4,  -40.48,  -25.68,  212.20,  -38.4];
% a_{bg} = [-51.22,  1.99, -1.009, -57.21,  -1.484,  -1.583,  -1.501, -267.3];

% Na23
% B_loc_target = [737.69];
% channel_series = ["aa"];
% ll_series = [0];
% B_loc = [738.223];
% Delta = [-159.8];
% a_{bg} = [-0.6757];
