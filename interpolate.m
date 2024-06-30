function y=interpolate(x,xk,yk)
Nx=max(size(xk));
rst=0;
for jj=1:Nx
    ljx=1;
    for ii=1:Nx
        if(abs(ii-jj)>0.001)
            ljx=ljx*(x-xk(ii))/(xk(jj)-xk(ii));
        end
    end
    rst=yk(jj)*ljx+rst;
end
y=rst;
end