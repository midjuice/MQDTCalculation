function y= ZeroEfun(Es,element,channel,B,ll,gJ,rs)

[Kc,~] =  Kcfunc(Es, element,channel,B,ll,gJ);
[~, ~, ~, ~, ~, ~, ~, ~, ~, mu, C6, ~] = Get_Params(element, channel);
beta6 = (C6*2*mu)^0.25;

nu0 =(2*ll+1)/4;
% from DOI: 10.1140/epjd/e2004-00127-x
x = (1/2)^(1/2) * (rs)^(-2);
fzero = (1/2)^(1/2) * (rs)^(1/2) * (cbesselj(nu0,x)*cos(pi*nu0/2) - cbessely(nu0,x)*sin(pi*nu0/2));
gzero = -(8)^(1/2) * (rs)^(1/2) * (cbesselj(nu0,x)*sin(pi*nu0/2) + cbessely(nu0,x)*cos(pi*nu0/2));


% From DOI: 10.1103/PhysRevA.64.010701
% x = (rs/beta6)^(-2)/2;
% fzero = 1/(2^(3/2)*cos(pi*nu0/2)) * rs^(1/2) * (cbesselj(nu0,x) + cbesselj(-nu0,x));
% gzero = -1/(2^(3/2)*sin(pi*nu0/2)) * rs^(1/2) * (cbesselj(nu0,x) - cbesselj(-nu0,x));


frow = fzero * ones(1,size(Kc,1));
grow = gzero * ones(1,size(Kc,1));

y = frow' + Kc * grow';

end



