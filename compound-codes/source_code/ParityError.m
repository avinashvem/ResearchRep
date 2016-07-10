function [K]=ParityError(UD,MDEL,Vcon,Ccon1)
E=zeros(1,max(Vcon(:)));
E(Vcon(MDEL,:))=1;
C1=sum(E(Ccon1),2)==size(Ccon1,2);

E=zeros(1,max(Vcon(:)));
E(Vcon(MDEL,:))=repmat(UD(MDEL),1,size(Vcon,2));
C2=mod(sum(E(Ccon1),2),2);
K=C1 & C2;