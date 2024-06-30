function y=cgfuncST(beta1,beta2,J,mj)
rst=0;
if(beta1==1 && beta2==1)
    if(J==1 && mj==1)
        rst=1;
    end
elseif(beta1==-1 && beta2==-1)
    if(J==1 && mj==-1)
        rst=1;
    end
elseif(beta1==1 && beta2==-1)
    if(mj==0)
        rst=1/sqrt(2);
    end
elseif(beta1==-1 && beta2==1)
    if(mj==0)
        if(J==1)
            rst=1/sqrt(2);
        else
            rst=-1/sqrt(2);
        end
    end
end
y=rst;
end