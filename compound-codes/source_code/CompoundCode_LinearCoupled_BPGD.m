clear; clc;
% For the LDPC part
dv1=3; % Variable-node degree
dc1=6; % Check-node degree

% For the LDGM Part
dv2=6; % Variable-node degree
dc2=3; % Check-node degree
dv=dv1+dv2;

c=1e2;
ms=dv2*(c); % Block Length / # of check-nodes in LDGM
ns=dc2*(c); % # of variable-nodes
m1s=ns*(dv1/dc1); % # of check-nodes in LDPC
e1s=ns*dv1; e2s=ns*dv2; es=ns*dv;

L=15; % Chain length of the SC system
w=3; % coupling parameter
Lw=L+w-1; % code-bits system length
ldpc_start=w; ldpc_end=L+w-1; 
ldpc_len=ldpc_end-ldpc_start+1; % ldpc-check system length
m=ms*L;
n=ns*Lw;
m1=m1s*ldpc_len;

rate=(dc2/dv2)*(1-(dv1/dc1));
dist_sha=finverse(@(x)(BinaryEntropy(x)),1-rate,0,0.5);
fprintf('rate=%f, Shannon distortion is %f\n',rate,dist_sha);

e1=Lw*e1s; e2=Lw*e2s; e=Lw*es;

E=zeros(1,e+1); % This is the main object in the sim

Vcon=zeros(n,dv);
Vcon(:,1:dv1)=reshape(1:e1,dv1,[])';
Vcon(:,dv1+1:dv)=e1+reshape(1:e2,dv2,[])';

P1=cell(Lw); P2=cell(Lw);

Ccon1=zeros(m1,dc1); Ccon2=zeros(m,dc2);

for i=1:Lw
    P1{i}=e1s*(i-1)+reshape(randperm(e1s),w,[]);
    P2{i}=e1+e2s*(i-1)+reshape(randperm(e2s),w,[]);
end

for i=ldpc_start:ldpc_end
    Etmp1=zeros(1,e1s);
    for j=1:w
       Ptmp1=P1{i-j+1};
       Etmp1((j-1)*e1s/w+1:j*e1s/w)=Ptmp1(w+1-j,:);
    end
    Etmp1=Etmp1(randperm(e1s));
    Ccon1((i-ldpc_start)*m1s+1:(i-ldpc_start+1)*m1s,:)=reshape(Etmp1,[],dc1);
end

for i=1:L
    Etmp2=zeros(1,e2s);
    for j=1:w
       Ptmp2=P2{i+j-1};
       Etmp2((j-1)*e2s/w+1:j*e2s/w)=Ptmp2(w+1-j,:);
    end
    Etmp2=Etmp2(randperm(e2s));
    Ccon2((i-1)*ms+1:i*ms,:)=reshape(Etmp2,[],dc2);
end

isFix4Cycles=0;
while(~isFix4Cycles)
    isFix2Cycles=0;
    while(~isFix2Cycles)
        [Vcon Ccon1 Ccon2 isFix2Cycles]=fix2cycles(Vcon,Ccon1,Ccon2,dv1,dv2);
    end
    [Vcon Ccon1 Ccon2 isFix4Cycles]=fix4cycles(Vcon,Ccon1,Ccon2,dv1,dv2);
end

%%%% Delete unused Vcon sockets
E_code=zeros(1,e+1);
E_code(Ccon1)=ones(size(Ccon1)); E_code(Ccon2)=ones(size(Ccon2));
Vcon(E_code(Vcon)==0)=e+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


M=1; % # of source sequences to average the distortion

beta=0.68; % inverse temperature
T=10; % Number of BP iterations before a forced decimation

DistProfile=zeros(L,M);
Dist=zeros(1,M);
tic;
for j=1:M
    UD=zeros(n,1);
    % Bernoulli symmetric source
    y=zeros(m,1); y(rand(m,1)<0.5)=1; % current instance of source seq.
    yr=y;
    E=zeros(1,e+1);
    count_n=n;
    MDEL=[]; % Variable bits deleted
    while(count_n > 0)
        for i=1:T
            % Variabel-node update
            E(end)=0; E(Vcon(MDEL,:))=0;
            E(Vcon)=repmat(sum(E(Vcon),2),1,dv1+dv2)-E(Vcon);
            E(E>20)=20; E(E<-20)=-20;

            % Check-node LDGM update
            E(Vcon(MDEL,dv1+1:dv))=Inf; 
            E(Ccon2)=(repmat(((-1).^yr)*tanh(beta),1,dc2)).*GenProd(E(Ccon2),1);
            E(Ccon2)=atanh(E(Ccon2));

            % Check-node LDPC update
            E(Vcon(MDEL,1:dv1))=repmat(((-1).^UD(MDEL))*Inf,1,dv1);
            E(Ccon1)=atanh(GenProd(E(Ccon1),1));
            E(E>20)=20; E(E<-20)=-20;
        end
        E(Vcon(MDEL,:))=0;
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
            parities_error_sequence=ParityError(UD,MDEL,Vcon,Ccon1);
            fprintf('%d/%d,err/tot=%d/%d\n',n-count_n,n,sum(parities_error_sequence),m1);
            parities_error_section=sum(reshape(parities_error_sequence,[],ldpc_len),1);
%             plot(sum(reshape(yr,[],L),1)/ms); ylim([0 0.6]); pause(0.1);
            if sum(parities_error_sequence)>0
                fprintf('Not Encoded Correctly!, parities in error=%d\n',sum(parities_error_sequence));
                disp(parities_error_section);
                fprintf('Reencoding...\n');
                yr=y;
                E=zeros(1,e+1);
                count_n=n;
                MDEL=[];
                UD=zeros(n,1);
            end
            toc;
        end
        if(count_n==0)
            E_code=zeros(1,e+1);
            E_code(Vcon)=repmat(UD,1,dv);
            parities_error_sequence=mod(sum(E_code(Ccon1),2),2);
            parities_error_section=sum(reshape(parities_error_sequence,[],ldpc_len),1);
%             plot(sum(reshape(yr,[],L),1)/ms); ylim([0 0.6]); pause(0.1);
            if sum(parities_error_sequence)==0
                fprintf('Encoding successful\n');
            else
                fprintf('Not Encoded Correctly!, parities in error=%d\n',sum(parities_error_sequence));
                disp(parities_error_section);
                fprintf('Reencoding...\n');
                yr=y;
                E=zeros(1,e+1);
                count_n=n;
                MDEL=[];
            end
        end
    end
    DistProfile(:,j)=sum(reshape(yr,[],L),1)/ms;
    DistDes=sort(DistProfile(:,j),'descend');
    Dist(j)=mean(DistDes(1:floor(end/2)));
    fprintf('Distortion for the curr. seq is %f, avg is %f \n',Dist(j), mean(Dist(1:j)) );
 end
