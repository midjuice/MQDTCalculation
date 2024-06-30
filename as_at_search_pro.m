clear
clc

element = 'Li6';
channel_series = ["ab", "ab", "ac", "ac", "aa", "ab", "bb", "bc"];
ll_series = [0, 0, 1, 0, 1, 1, 1, 0];
left_series = [600, 520, .1, 500, 100, 100, 200, 700];
right_series = [1000, 570, 400, 1000, 300, 300, 400, 1000];
points_series = [500, 1000, 500, 500, 500, 500, 500, 500];
B_loc_target = [834.1, 543, 225, 690.4, 159, 185, 215, 811];
gJ = 2.0023;


% scattering lengths at the beginning
KcST_start = [4, 20]; %[14.4119, 0.3981]; %[17.7924147429038, 0.396695033941700];

% numbers of iteration
num_iter = 5;

% standard deviation of variation of the two parameters
sigma = [0.01, 0.01];

% tolerance to end the for loop
tol = 1;

% initialization
KcST = KcST_start;
Ks = KcST(1);
Kt = KcST(2);
B_loc = B_loc_func(element, channel_series, ll_series, left_series, right_series, points_series, gJ, Ks, Kt);
samples = zeros(num_iter, length(KcST)+length(channel_series)+1);

% Optimization loop
for i = 1:num_iter
    KcST_new = KcST + randn(1, length(KcST)).*sigma;
    Ks_new = KcST_new(1);
    Kt_new = KcST_new(2);
    B_loc_new = B_loc_func(element, channel_series, ll_series, left_series, right_series, points_series, gJ, Ks_new, Kt_new);
    % accept rate
    % alpha = exp(-abs(y_new - y_target));
    dist = pdist([B_loc_target; B_loc]);
    dist_new = pdist([B_loc_target; B_loc_new]);
    if dist_new < dist
        KcST = KcST_new;
        B_loc = B_loc_new;
    end
    if dist_new < tol
        break;
    end   
    % save the parameters
    samples(i, 1:2) = KcST;
    samples(i, 3:end-1) = B_loc;
    samples(i, end) = dist_new;
    % report the progress
    disp(i);
end


function y = B_loc_func(element, channel_series, ll_series, left_series, right_series, points_series, gJ, Ks, Kt)
    num_resonance = length(channel_series);
    y = zeros(1, num_resonance);
    parfor s=1:num_resonance
        channel = channel_series(s);
        ll = ll_series(s);
        left = left_series(s);
        right = right_series(s);
        points = points_series(s);

        [Ia, mf1k, alpha1k, mf2k, alpha2k, gI, sgn, Ehf, E6, ~] = Get_Params(element, channel);
        Bs = linspace(left, right, points);
        output = zeros(1, points);
        
        for i=1:points
            B = Bs(i);
            output(i) = real((-1)^ll+1/KeffuncID(0,ll,B,mf1k,alpha1k,mf2k,alpha2k,gJ,gI,Ehf,Ia,Ks,Kt,E6,sgn));
        end
        [~, Index] = max(output);
        y(s) = Bs(Index);
    end
end

