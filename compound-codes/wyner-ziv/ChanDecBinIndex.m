function [UD] = ChanDecBinIndex(Vcon,Ccon1,Ccon2,BinIndex,err,y,L,w)

y=1-2*y;
n=size(Vcon,1); ns=n/(L+w-1); dv=size(Vcon,2);
m1=size(Ccon1,1); dc1=size(Ccon1,2);
dc2=size(Ccon2,2);
llr0=log((1-err)/err);
iter=300;
e=max(max(max(Vcon(:)),max(Ccon1(:))),max(Ccon2(:)));
assert(numel(BinIndex)==m1);

E=zeros(1,e); Eold=zeros(size(E));
for i=0:iter
    %Variable-node update
    E(end)=0;
    E(Vcon)=repmat(sum(E(Vcon),2),1,dv)-E(Vcon);
    E(Vcon(1:((w-1)*ns),:))=20; % boundary condition in a SC system
    E(E>20)=20; E(E<-20)=-20;
    E(end)=0;

    % Check-node update
    E(end)=Inf;
    E(Ccon1)=repmat((-1).^BinIndex,1,dc1).*GenProd(E(Ccon1),1/2);
    E(Ccon2)=repmat(tanh(llr0/2).*y,1,dc2).*GenProd(E(Ccon2),1/2);
    E=2*atanh(E);
    E(E>20)=20; E(E<-20)=-20;
    
    if(max(abs(Eold(1:e)-E(1:e)))<1e-6)
        fprintf('Messages Converged! Stopping...\n');
        break;
    end
    Eold=E;
end
E(end)=0;
UD=sum(E(Vcon),2)<0; UD(1:(w-1)*ns)=0;

return;
