clear
clc

element = 'Li6';
channel_series = ["ab", "ab", "ac", "ac", "aa", "ab", "bb", "bc"];
ll_series = [0, 0, 1, 0, 1, 1, 1, 0];
left_series = [600, 520, .1, 500, 100, 100, 200, 700];
right_series = [1000, 570, 400, 1000, 300, 300, 400, 1000];
points_series = [500, 1000, 500, 500, 500, 500, 500, 500];
B_loc_target = [834.1, 543, 225, 690.4, 159, 185, 215, 811];
B_loc_width = [300, 0.1, 10, 122, 159, 28, 42, 200];
weight = (1./B_loc_width)./sum(1./B_loc_width);
gJ = 2.0023;

% scattering lengths at the beginning
a0_start = [44.7634 -2134.8]; %[44.36364399, -1995.706848];

% numbers of iteration
num_iter = 30;

% standard deviation of variation of the two parameters
sigma = [0.1, 100];

% tolerance to end the for loop
tol = 1;

% initialization
a0 = a0_start;
B_loc = B_loc_func(element, channel_series, ll_series, left_series, right_series, points_series, gJ, a0);
samples = zeros(num_iter, length(a0)+length(channel_series)+1);

% Optimization loop
for i = 1:num_iter
    a0_new = a0 + randn(1, length(a0)).*sigma;
    B_loc_new = B_loc_func(element, channel_series, ll_series, left_series, right_series, points_series, gJ, a0_new);
    dist = dist_func(B_loc_target, B_loc, weight);
    dist_new = dist_func(B_loc_target, B_loc_new, weight);
    if dist_new < dist
        a0 = a0_new;
        B_loc = B_loc_new;
    end
    if dist_new < tol
        break;
    end   
    % save the parameters
    samples(i, 1:2) = a0;
    samples(i, 3:end-1) = B_loc;
    samples(i, end) = dist;
    % report the progress
    disp(i);
end


function y = B_loc_func(element, channel_series, ll_series, left_series, right_series, points_series, gJ, a0)
    num_resonance = length(channel_series);
    y = zeros(1, num_resonance);
    parfor s=1:num_resonance
        channel = channel_series(s);
        ll = ll_series(s);
        left = left_series(s);
        right = right_series(s);
        points = points_series(s);

        [Ia, mf1k, alpha1k, mf2k, alpha2k, gI, sgn, Ehf, E6, beta6, ~] = Get_Params(element, channel);
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
        y(s) = Bs(Index);
    end
end

function dist = dist_func(B_loc_target, B_loc, weight)
    dist = sum((B_loc_target-B_loc).^2.*weight);
end

