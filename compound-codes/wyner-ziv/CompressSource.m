function [Dist DistProfile UD yr]=CompressSource(Vcon,Ccon,y)

global L w;

n=size(Vcon,1); % # of code-bits
m=size(Ccon,1); % # of check-nodes
dv=size(Vcon,2); % code-bit degree
dc=size(Ccon,2); % check-node degree
ns=n/(L+w-1);
beta=1.04; % inverse temperature
T=10; % Number of BP iterations before a forced decimation

e=max(max(Vcon(:)),max(Ccon(:)));
UD=zeros(n,1);
yr=y;
E=zeros(1,e);
count_n=n;
MDEL=[]; % Variable bits deleted

% Puncturing code-bits in left-most w-1 sections
MDEL=[MDEL 1:(ns*(w-1))];
UD(1:(ns*(w-1)))=0;
count_n=count_n-(ns*(w-1));

while(count_n > 0)
    for i=1:T
        % Variabel-node update
        E(end)=0;
        E(Vcon)=repmat(sum(E(Vcon),2),1,dv)-E(Vcon);

        % Check-node update
        E(end)=Inf; E(Vcon(MDEL,:))=Inf; 
        E(Ccon)=(repmat(((-1).^yr)*tanh(beta),1,dc)).*GenProd(E(Ccon),1);
        E(Ccon)=atanh(E(Ccon));
        E(E>20)=20; E(E<-20)=-20;

        E(Vcon(MDEL,:))=0; E(end)=0;
    end
    % Possibility of decimation?     
    ND=1:n; ND(MDEL)=[];
    idlow=ceil(ND(1)/ns); idhigh=idlow+w-1;
    idlow=(idlow-1)*ns+1; idhigh=idhigh*ns;
    ND=ND(ND>=idlow & ND<=idhigh);

    bt=sum(E(Vcon),2);
    [B,i]=max(abs(bt(ND)));
    nd=ND(i);

    if(B > 0)
%         ud=rand>(1+tanh(bt(nd)))/2;
        ud=bt(nd)<0;
        assert(all(MDEL-nd));
    else
        ND=1:n; ND(MDEL)=[];
        nd=ND(1);
        ud=(rand<0.5);
    end
    UD(nd)=ud;
    MDEL=[MDEL nd];
    count_n=count_n-1;
    if(ud==1)
        E_code=zeros(e,1);
        E_code(Vcon(nd,:))=1;
        yr=mod(sum([E_code(Ccon) yr],2),2);
    end
    if(mod(count_n,100)==0)
        fprintf('Finished %d/%d code-bits\n',n-count_n,n);
    end
end

DistProfile=sum(reshape(yr,[],L),1)/(m/L);
Dist=mean(DistProfile(w+2:L-w));
