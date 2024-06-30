clear
clc

element = 'Li6';
channel = 'ac';
ll = 1;
gJ = 2.0023;

left = .1;
right = 600;
points = 1000;


% 设置初始值
a0_start = [44, -2235];

% 设置迭代次数
num_iter = 20;

% 设置标准差
sigma = [0.1, 20];

% 设置偏差容忍度
tol = 0.1;

% 初始化参数
a0 = a0_start;
y = B_loc_func(element, channel, ll, gJ, left, right, points, a0);
y_target = 225;
samples = zeros(num_iter, 3);


% 迭代 MCMC
for i = 1:num_iter
    % 生成新的参数
    a0_new = a0 + randn(1, length(a0)).*sigma;
    
    % 计算新的函数值
    y_new = B_loc_func(element, channel, ll, gJ, left, right, points, a0_new);
    % 计算接受率
    % alpha = exp(-abs(y_new - y_target));
    
    % 判断是否接受新的参数
    if abs(y_new-y_target) < abs(y-y_target)
        a0 = a0_new;
        y = y_new;
    end
    
    % 判断是否满足偏差容忍度
    if abs(y - y_target) < tol
        break;
    end
    
    % 保存参数
    samples(i, 1:2) = a0;
    samples(i, 3) = y;
    disp(i);
end


function y = B_loc_func(element, channel, ll, gJ, left, right, points, a0)
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
    y = (Bs(Index)+Bs(Index+1))/2;
end

