function [x]=finverse(f,y,a,b)

eps=1e-6;
if(f(b-eps) > f(a+eps))
    xlow=a; xhigh=b;
    while(xhigh-xlow > eps)
        xmid=(xlow+xhigh)/2;
        if(y>f(xmid))
            xlow=xmid;
        elseif(y<f(xmid))
            xhigh=xmid;
        else
            x=xmid;
            return;
        end
    end
elseif(f(b-eps) < f(a+eps))
    xlow=b; xhigh=a;
    while(xhigh-xlow > eps)
        xmid=(xlow+xhigh)/2;
        if(y>f(xmid))
            xhigh=xmid;
        elseif(y<f(xmid))
            xlow=xmid;
        else
            x=xmid;
            return;
        end
    end
end
x=xmid;