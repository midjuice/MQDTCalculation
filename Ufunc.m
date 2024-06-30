% alpha denotes (M_F)_{1 or 0 or -1}
% beta denotes M_j1 or M_j2 (1/2 or -1/2)
% mf values -> -I_A-1/2, -I_A, I_A+1/2, ... , I_A, I_A+1/2

function y=Ufunc(alpha,beta,mf,Ia,x)
rst=0;
if(2*mf==2*Ia+1)
    if(alpha==0 && beta==1)
        rst=1;
    end
elseif(-2*mf==2*Ia+1)
    if(alpha==0 && beta==-1)
        rst=1;
    end
else
    temp=2*mf/(2*Ia+1);
    theta=atan(sqrt(1-temp^2)/(temp+x));
    if(theta<0)
        theta=theta+pi;
    end
    if(alpha*beta>0)
        rst=sin(theta/2);
    else
        if(alpha>0)
            rst=-cos(theta/2);
        else
            rst=cos(theta/2);
        end
    end
end
y=rst;
end