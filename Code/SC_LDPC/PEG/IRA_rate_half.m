function [Vcon,Ccon,Emax]=IRA_rate_half(n,repeat_flag)
% n: Number of bit nodes
global repeat
% n=200;
k=n/2;
dv=3;
dc=5;

Emax=k*dv+2*k-1;
Vcon=(Emax+1)*ones(n,dv);
Ccon=(Emax+1)*ones(k,dc);

Vcon(1:k,1:dv)=reshape(1:k*dv,[],dv);
Ccon(1:k,1:dv)=reshape(randperm(k*dv),[],dv);
e=k*dv;

Vcon(k+1,1)=e+1; Ccon(1,dv+1)=e+1;
e=e+1;
for i=2:k
    Vcon(k+i,1:2)=e+1:e+2;
    Ccon(i,dv+1:dv+2)=e+1:e+2;
    e=e+2;
end
Emax=e;

    
if repeat_flag
    N=size(Vcon,1); m=size(Ccon,1);
    Emax_old=Emax;
    Emax=Emax+2*(repeat-1)*N;
    Vcon_old=Vcon;
    
    Vcon=(Emax+1)*ones(N*repeat,size(Vcon_old,2)+1);
    Vcon(1:N,1:size(Vcon_old,2))=Vcon_old;
        
    Vcon(Vcon==Emax_old+1)=Emax+1;
    Ccon(Ccon==Emax_old+1)=Emax+1;

    e=Emax_old;
    for i=1:N
        Vcon(i,end)=e+1;
        for j=1:repeat-2
            Vcon(i+j*N,1:2)=e+2:e+3;
        end
        Vcon(i+(repeat-1)*N,1)=e+2*(repeat-2)+1+1;
        Ccon(m+1:m+repeat-1,1:2)=reshape(e+1:e+2*(repeat-1),[],2);
        Ccon(m+1:m+repeat-1,3:end)=Emax+1;
        
        m=m+repeat-1;
        e=e+2*(repeat-1);
    end 
end
isFix2Cycles=1;
while isFix2Cycles>0
    [Vcon Ccon isFix2Cycles]=fix2cycles(Vcon,Ccon);
end

isFix4Cycles=1;
while isFix4Cycles>0
    [Vcon Ccon isFix4Cycles]=fix4cycles(Vcon,Ccon);
end

end
