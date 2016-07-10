function [K]=ParityError(UD,MDEL,Vcon,Ccon1,BinIndex)
e=max(max(Vcon(:)),max(Ccon1(:)));
E=zeros(1,e);
E(Vcon(MDEL,:))=1;
C1=sum(E(Ccon1),2)==size(Ccon1,2);

E=zeros(1,e);
E(Vcon(MDEL,:))=repmat(UD(MDEL),1,size(Vcon,2));
C2=mod(mod(sum(E(Ccon1),2),2)+BinIndex,2);
K=C1 & C2;
