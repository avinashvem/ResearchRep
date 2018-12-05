function Lout=LLR_real_adder(L,y,sig,G)
% sig=AWGN noise variance
% y- Nx1 dim. signal rx from U-user BPSK input real adder GMAC channel
% L- NxU llrs of the connected bit nodes

llr_max=20;
[N,U]=size(L);
Lout=zeros(size(L));
U0=(numel(G(1,:))-1)/2;
% G(i,1) corresponds to Q(y(i),-U0)
% G(i,j) corresponds to Q(y(i),-U0+j-1)
for i=1:N
    idx=find(imag(L(i,:))==0);
    Li=L(i,idx);
    Ui=numel(Li);
    
    if Ui==0
       Lout(i,1:U)=2*y(i)/sig^2;      
    
    elseif Ui==1
       l1=exp(Li(1));
       Lout(i,1:U)=log((1+exp(l1+2*(y(i)-1)/sig^2))./(exp(l1)+exp(-2*(y(i)+1)/sig^2)));
       Lout(i,idx(1))=2*y(i)/sig^2;      
  
    elseif Ui==2
       Li=exp(Li);
       sumD=[1 sum(Li) prod(Li)];
       Lout(i,1:U)=log(sum(G(i,U0-Ui+2:2:U0+Ui+2).*sumD)/sum(G(i,U0-Ui:2:U0+Ui).*sumD)); 
       
       l1=Li(1); l2=Li(2);
       Lout(i,idx(1))= log(sum([1 l2].*G(i,U0+1:2:U0+3))/sum([1 l2].*G(i,U0-1:2:U0+1)));
       Lout(i,idx(2))= log(sum([1 l1].*G(i,U0+1:2:U0+3))/sum([1 l1].*G(i,U0-1:2:U0+1)));

    elseif Ui==3
       Li=exp(Li);
       l1=Li(1); l2=Li(2); l3=Li(3);
       sumD=[1 sum(Li) l1*l2+l1*l3+l2*l3 prod(Li)];
       Lout(i,1:U)=log(sum(G(i,U0-Ui+2:2:U0+Ui+2).*sumD)/sum(G(i,U0-Ui:2:U0+Ui).*sumD)); 
       
       Lout(i,idx(1))= log(sum([1 l2+l3 l3*l2].*G(i,U0:2:U0+4))/sum([1 l2+l3 l3*l2].*G(i,U0-2:2:U0+2)));
       Lout(i,idx(2))= log(sum([1 l1+l3 l1*l3].*G(i,U0:2:U0+4))/sum([1 l1+l3 l1*l3].*G(i,U0-2:2:U0+2)));
       Lout(i,idx(3))= log(sum([1 l1+l2 l1*l2].*G(i,U0:2:U0+4))/sum([1 l1+l2 l1*l2].*G(i,U0-2:2:U0+2)));
    
    else
        
       Li=exp(Li);
       sumD=[1 Li(1)];
       for u=2:Ui
         sumD=conv(sumD,[1 Li(u)]);
       end
       Lout(i,1:U)=log(sum(G(i,U0-Ui+2:2:U0+Ui+2).*sumD)/sum(G(i,U0-Ui:2:U0+Ui).*sumD));
        
       for u=1:Ui
         [extrD,r]=deconv(sumD,[1 Li(u)]);
         extrD=max(0,extrD);
         Lout(i,idx(u))=log(sum(G(i,U0-Ui+3:2:U0+Ui+1).*extrD)/sum(G(i,U0-Ui+1:2:U0+Ui-1).*extrD));
       end
    end
   Lout(i,isnan(Lout(i,:)))=0;
end
Lout(Lout>llr_max)=llr_max;
Lout(Lout<-llr_max)=-llr_max;
Lout(isnan(Lout))=0;
end