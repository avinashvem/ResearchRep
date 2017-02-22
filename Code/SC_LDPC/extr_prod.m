function Eop=extr_prod(E,Ccon,Parity)
dc=length(Ccon(1,:));

E=tanh(E/2);
Eop=E;

chkProd=prod([E(Ccon) tanh(Parity/2)],2);

if prod(chkProd)~=0
    Eop(Ccon)=repmat(chkProd,1,dc)./E(Ccon);
else
    nonZeroIdx= chkProd~=0;
    if ~isempty(find(nonZeroIdx,1))
    Eop(Ccon(nonZeroIdx,:))=repmat(chkProd(nonZeroIdx),1,dc)./ E(Ccon(nonZeroIdx,:));
    end
    zeroIdx=find(nonZeroIdx==0);
    for i=1:length(zeroIdx)
        cumProd=+1*sign(Parity(zeroIdx(i)));
        zeroEdges=[];
        for j=1:dc
            if E(Ccon(zeroIdx(i),j))~=0
             cumProd=cumProd*E(Ccon(zeroIdx(i),j)); 
            else
                zeroEdges=[zeroEdges j];
            end
        end
        
        if length(zeroEdges)>1 %Equiv to More than 1 erasure
          Eop(Ccon(zeroIdx(i),:))= 0;
        elseif ~isempty(zeroEdges) 
            Eop(Ccon(zeroIdx(i),:))= 0;
            Eop(Ccon(zeroIdx(i),zeroEdges(1)))= cumProd;
        else
            error('Throw error');
        end
    end
end