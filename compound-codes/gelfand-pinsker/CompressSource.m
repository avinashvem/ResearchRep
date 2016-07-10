function [Dist DistProfile UD yr]=CompressSource(Vcon,Ccon1,Ccon2,y,MSG_I,M)

% CompressSource(Vcon,Ccon1,Ccon2,S,MSG_I,M)

global L w;

n=size(Vcon,1); % # of code-bits
m1=size(Ccon1,1); % # of LDPC check-nodes
m=size(Ccon2,1); % # of LDGM check-nodes
dv=size(Vcon,2); % code-bit degree
dc1=size(Ccon1,2); % LDPC check-node degree
dc2=size(Ccon2,2); % LDGM check-node degree
ldpc_len=L-w+1;
ns=n/(L+w-1);
beta=0.7; % inverse temperature
T=10; % Number of BP iterations before a forced decimation

e=max(max(max(Vcon(:)),max(Ccon1(:))),max(Ccon2(:)));
UD=zeros(n,1);
yr=y;
E=zeros(1,e);
count_n=n;
MDEL=[]; % Variable bits deleted

% Puncturing code-bits in left-most w-1 sections
MDEL=[MDEL 1:(ns*(w-1))];
UD(1:(ns*(w-1)))=0;
count_n=count_n-(ns*(w-1));

%Effective Bin Index in LDPC Checks
BinIndex=zeros(m1,1);
BinIndex(MSG_I)=M;

while(count_n > 0)
    for i=1:T
        % Variabel-node update
        E(end)=0;
        E(Vcon)=repmat(sum(E(Vcon),2),1,dv)-E(Vcon);
        E(E>20)=20; E(E<-20)=-20;

        % Check-node LDGM update
        E(Vcon(MDEL,:))=Inf; 
        E(Ccon2)=(repmat(((-1).^yr)*tanh(beta),1,dc2)).*GenProd(E(Ccon2),1);
        E(Ccon2)=atanh(E(Ccon2));

        % Check-node LDPC update
        E(Vcon(MDEL,:))=repmat(((-1).^UD(MDEL))*Inf,1,dv);
        E(Ccon1)=repmat((-1).^BinIndex,1,dc1).*atanh(GenProd(E(Ccon1),1));
        E(E>20)=20; E(E<-20)=-20;

        E(Vcon(MDEL,:))=0;
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
        ud=rand>(1+tanh(bt(nd)))/2;
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
        E_code=zeros(e+1,1);
        E_code(Vcon(nd,:))=1;
        yr=mod(sum([E_code(Ccon2) yr],2),2);
    end
    if(mod(count_n,100)==0)
        parities_error_sequence=ParityError(UD,MDEL,Vcon,Ccon1,BinIndex);
        fprintf('%d/%d,err/tot=%d/%d\n',n-count_n,n,sum(parities_error_sequence),m1);
        parities_error_section=sum(reshape(parities_error_sequence,[],ldpc_len),1);
        if sum(parities_error_sequence)>0
            fprintf('Not Encoded Correctly!, parities in error=%d\n',sum(parities_error_sequence));
            disp(parities_error_section);
            fprintf('Reencoding...\n');
            yr=y;
            E=zeros(1,e);
            MDEL=1:(ns*(w-1));
            UD=zeros(n,1);
            UD(1:(ns*(w-1)))=0; % puncturing
            count_n=n-(ns*(w-1));
        end
    end
    if(count_n==0)
        E_code=zeros(1,e);
        E_code(Vcon)=repmat(UD,1,dv);
        parities_in_error=sum(mod(sum(E_code(Ccon1),2)+BinIndex,2));
        assert(parities_in_error==0);
        fprintf('Encoding successful\n');
    end
end
DistProfile=sum(reshape(yr,[],L),1)/(m/L);
Dist=mean(DistProfile(w+2:L-w));

