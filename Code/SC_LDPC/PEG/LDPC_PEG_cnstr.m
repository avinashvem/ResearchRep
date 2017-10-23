function [Vcon,Ccon,Emax]=LDPC_PEG_cnstr(n,m,L,R,repeat_flag)
% n: Number of bit nodes
% m: # check nodes
% L(1,:): left d.d L(x)=\sum L(1,i) x^ldeg(i)
%ldeg=L(2,:)

global repeat
    ldeg=L(2,:);  rdeg=R(2,:);
    L=L(1,:); L=round(n*L);
    R=R(1,:); R=round(m*R);
    
    n=sum(L); m=sum(R);
    cnt=0;
    for i=1:length(L)
      ndeg_vec(cnt+1:cnt+L(i))=ldeg(i);
      cnt=cnt+L(i);
    end
    cnt=0;
    for i=1:length(R)
      mdeg_vec(cnt+1:cnt+R(i))=rdeg(i);
      cnt=cnt+R(i);
    end
    %{
    [Vcon_C,Ccon_V,girth]=ldpc_peg_const(n,m,d);
    
    Emax=numel(find(Vcon_C>0));
    Vcon=zeros(size(Vcon_C));
    Vcon(Vcon_C>0)=1:Emax;
     
    d=zeros(n,1);
    Ccon=(Emax+1)*ones(m,size(Ccon_V,2));
    for i=1:m
        j=1;
        while j<=size(Ccon_V,2) && Ccon_V(i,j)>0
            Vij=Ccon_V(i,j);
            Ccon(i,j)=Vcon(Vij,d(Vij)+1);
            d(Vij)=d(Vij)+1;
            j=j+1;
        end
    end
   
    Vcon(Vcon==0)=Emax+1;
   %}
    assert(sum(ndeg_vec)==sum(mdeg_vec))
    [Vcon,Ccon,Emax]=ldpc_peg_constr_michael(n,m,ndeg_vec,mdeg_vec);
    
if repeat_flag
%     global repeat
    N=size(Vcon,1); m=size(Ccon,1);
    Emax_old=Emax;
    Emax=Emax+2*(repeat-1)*N;
    Vcon_old=Vcon;
%     Ccon_old=Ccon;
    
    Vcon=(Emax+1)*ones(N*repeat,size(Vcon_old,2)+1);
    Vcon(1:N,1:size(Vcon_old,2))=Vcon_old;
    
%     Ccon=(Emax+1)*ones(m*repeat,size(Vcon_old,2)+1);
%     Ccon(1:N,1:size(Vcon_old,2))=Vcon_old;
    
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
end