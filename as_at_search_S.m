clear
clc

element = 'Li6';
channel_series = ["ab", "ab", "ac", "ac", "aa", "ab", "bb", "bc"];
ll_series = [0, 0, 1, 0, 1, 1, 1, 0];
left_series = [600, 530, .1, 500, 100, 100, 200, 700];
right_series = [1000, 560, 400, 1000, 300, 300, 400, 1000];
points_series = [500, 500, 500, 500, 500, 500, 500, 500];
B_loc_target = [834.1, 543, 225, 690.4, 159, 185, 215, 811];
gJ = 2.0023;


% x denotes a series of Ks, while y denotes a series of Kt
% KcST_start = [17.7924147429038, 0.396695033941700];
x = 15:0.1:18;
y = 0.35:0.01:0.45;
[X, Y] = meshgrid(x, y);

% numbers of iteration
num_iterx = length(x);
num_itery = length(y);

% standard deviation of variation of the two parameters
% sigma = [0.002, 0.002];

% initialization
% B_loc = B_loc_func(element, channel_series, ll_series, left_series, right_series, points_series, gJ, KcST);
samples = zeros(num_itery, num_iterx, 11);

for i = 1:num_iterx
    for j = 1:num_itery
        Ks = X(j, i);
        Kt = Y(j, i);
        B_loc = B_loc_func(element, channel_series, ll_series, left_series, right_series, points_series, gJ, Ks, Kt);
        dist = pdist([B_loc_target; B_loc]);
        % save the parameters
        samples(j, i, 1:2) = [Ks, Kt];
        samples(j, i, 3) = dist;
        samples(j, i, 4:end) = B_loc;
    end
    disp(i);
end

surf(X, Y, samples(:,:,3));
xlabel("Ks");
ylabel("Kt");
zlabel("Distance");
[row, col] = find(samples(:,:,3)==min(min(samples(:,:,3))));
disp(samples(row, col, :))


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

