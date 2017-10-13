% function [Vcon,Ccon,Emax]=ECon_LDPC(regular_degree)
function [Vcon,Ccon,Emax]=ECon_LDPC(n,m,L,R,repeat_flag)
% For irregular, currently builds the graph only for
% LDPC code d.d. of Right-reg(can be left-irregular) works
     %{
if regular_degree

    if mod(M*dv/dc,1)==0
        Vtot=M;
        Ctot=M*dv/dc;
    else 
        G=gcd(dc,dv);
        m=dv/G; n=dc/G; % It should be m*dc=n*dv
        Vtot=M*n;  Ctot=M*m;
    end

    Emax=Vtot*dv;

    Vcon=reshape(1:Emax,dv,[])';
    Ccon=reshape(randperm(Emax),dc,[])';
else
%}     

    ldeg=L(2,:);  rdeg=R(2,:);
    L=L(1,:); L=round(n*L);
    R=R(1,:); R=round(m*R);
    
    n=sum(L); m=sum(R);
    %{
    cnt=0;
    for i=1:length(L)
      ndeg_vec(icnt+1:cnt+L(i))=ldeg(i);
      cnt=cnt+L(i);
    end
    cnt=0;
    for i=1:length(R)
      mdeg_vec(cnt+1:cnt+R(i))=rdeg(i);
      cnt=cnt+R(i);
    end
%}
    assert(sum(L.*ldeg)==sum(R.*rdeg))
    Emax=sum(L.*ldeg);
    %{
    L=[0.315/2 0.2422/3 0.1/17]; L=L./sum(L);
    ldeg=[2 3 17];
    R=1; rdeg=8;
    %}
    
    Vcon=(Emax+1)*ones(n,max(ldeg));
    e=0;bits=0;
    for i=1:length(ldeg)
       Vcon(bits+1:bits+L(i),1:ldeg(i))=reshape(e+1:e+L(i)*ldeg(i),[],ldeg(i));
       bits=bits+L(i);
       e=e+L(i)*ldeg(i);
    end
    
    Ccon=(Emax+1)*ones(m,max(rdeg));
    edges=randperm(Emax);
    e=0;chks=0;
    for i=1:length(rdeg)
       Ccon(chks+1:chks+R(i),1:rdeg(i))=reshape(edges(e+1:e+R(i)*rdeg(i)),[],rdeg(i));
       chks=chks+R(i);
       e=e+R(i)*rdeg(i);
    end
    
    
if repeat_flag
    global repeat
    N=size(Vcon,1); m=size(Ccon,1);
    Emax_old=Emax;
    Emax=Emax+2*(repeat-1)*N;
    Vcon_old=Vcon;
    
    Vcon=(Emax+1)*ones(N*repeat,size(Vcon_old,2)+1);
    Vcon(1:N,1:size(Vcon_old,2))=Vcon_old;
    
    Vcon(Vcon==Emax_old+1)=Emax+1;
    Ccon(Ccon==Emax_old+1)=Emax+1;
    Ccon=[Ccon;(Emax+1)*ones(N*(repeat-1),size(Ccon,2))];

    e=Emax_old;
    for i=1:N
        idx=find(Vcon(i,:)==Emax+1,1,'first');
        Vcon(i,idx)=e+1;
        for j=1:repeat-2
            Vcon(i+j*N,1:2)=e+2:e+3;
        end
        Vcon(i+(repeat-1)*N,1)=e+2*(repeat-2)+1+1;
        Ccon(m+1:m+repeat-1,1:2)=reshape(e+1:e+2*(repeat-1),[],2);
        
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