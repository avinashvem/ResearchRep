clear; clc;

% For the LDP part
dv1=3; % Variable-node degree
dc1=6; % Check-node degree

% For the LDGM Part
dv2=6; % Variable-node degree
dc2=3; % Check-node degree
dv=dv1+dv2;

c=2e2;
ms=dv2*(c); % Block Length / # of check-nodes in LDGM
ns=dc2*(c); % # of variable-nodes
m1s=ns*(dv1/dc1); % # of check-nodes in LDPC

L=10; % Chain length of the SC system
w=2; % coupling parameter

m=ms*L; % Total # of check-nodes in LDGM / Block length
n=ns*L; % Total # of variable-nodes
m1=m1s*L; % Total # of check-nodes in LDPC

e1s=dv1*ns; % # of edges in LDPC part in single system
e2s=dv2*ns; % # of edges in the LDGM part in single system
es=e1s+e2s; % Total # of edges in single system

e1=e1s*L;
e2=e2s*L;
e=es*L; % Total # of edges

rate=(dc2/dv2)*(1-(dv1/dc1)); 
% dist_sha=finverse(@(x)(BinaryEntropy(x)),1-rate,0,0.5);
% fprintf('rate=%f, Shannon distortion is %f\n',rate,dist_sha);
fprintf('Block length  = %d\n',m);

E=zeros(1,e); % This is the main object in the sim
E_check_code=zeros(1,e); % Object to check if the code-bits satisfy LDPC checks

Vcon=zeros(n,dv1+dv2);
for i=1:L
    Vcon((1+(i-1)*ns) : (i*ns),:)= (i-1)*es + [reshape(1:e1s,dv1,[])' reshape(e1s+1:e1s+e2s,dv2,[])'];
end

% Constructing a Spatially-Coupled Code
Ccon1=zeros(m1,dc1); Ccon2=zeros(m,dc2);
P1=cell(L); P2=cell(L);

for i=1:L
    P1{i}=es*(i-1)+reshape(randperm(e1s),w,[]); % w bundles of variable-node sockets in each system
    P2{i}=es*(i-1)+reshape(e1s+randperm(e2s),w,[]);
end
for i=1:L
    Etmp1=zeros(1,e1s); Etmp2=zeros(1,e2s);
    for j=1:w
        if(i+j-1>L)
            eff_ind=i+j-1-L;
        else
            eff_ind=i+j-1;
        end
        Ptmp1=P1{eff_ind}; Ptmp2=P2{eff_ind};
        Etmp1((1+(j-1)*e1s/w) : (j*e1s/w))=Ptmp1(w+1-j,:);
        Etmp2((1+(j-1)*e2s/w) : (j*e2s/w))=Ptmp2(w+1-j,:);
    end
    Etmp1=Etmp1(randperm(e1s)); Etmp2=Etmp2(randperm(e2s));
    
    Ccon1(((i-1)*m1s+1) : (i*m1s) , :)=reshape(Etmp1,[],dc1);
    Ccon2(((i-1)*ms+1) : (i*ms) , :)=reshape(Etmp2,[],dc2);
end
imid=ceil(L/2);
ilow0=imid-w/2+1; ihigh0=mod(ceil(imid+w/2),L);

M=1; % # of source sequences to average the distortion

beta=0.8; % inverse temperature
T=10; % Number of BP iterations before a forced decimation
alp=4.25; % threshold for decimation

DistProfile=zeros(L,M);
Distortion=zeros(1,M);

for j=1:M
    UD=zeros(n,1);
    % Bernoulli symmetric source
    y=zeros(m,1); y(rand(m,1)<0.5)=1; % current instance of source seq.
    yr=y;
    E=zeros(1,e);
    count_n=n;
    MDEL=[]; % Variable bits deleted
    prev_parities_error=0;
    while(count_n > 0)
        for i=1:T
            % Variabel-node update
            E(Vcon)=repmat(sum(E(Vcon),2),1,dv1+dv2)-E(Vcon);

            % Check-node LDGM update
            E(Vcon(MDEL,dv1+1:dv))=Inf; 
            E(Ccon2)=(repmat(((-1).^yr)*tanh(beta),1,dc2)).*GenProd(E(Ccon2),beta);
            E(Ccon2)=(1/beta)*atanh(E(Ccon2));
            
            % Check-node LDPC update
            E(Vcon(MDEL,1:dv1))=repmat(((-1).^UD(MDEL))*Inf,1,dv1);
            E(Ccon1)=(1/beta)*atanh(GenProd(E(Ccon1),beta));
            E(E>20)=20; E(E<-20)=-20;

            E(Vcon(MDEL,:))=0;
            
            % Possibility of decimation?
            bt=beta*sum(E(Vcon),2);
            B=max(abs(bt));
            if(B>alp)
               break;
            end
        end
        if(B > 0)
            RAND_DECIM=0;
            I=find(abs(bt)==B);
            nd=I(randi(length(I),1,1));
            ud=bt(nd)<0;
            assert(all(MDEL-nd));
        else
            RAND_DECIM=1;
            idlow=(ilow0-1)*ns+1; idhigh=ihigh0*ns;
            MDEL1=MDEL(MDEL>=idlow & MDEL<=idhigh)-idlow+1;
            ND=idlow:idhigh;
            if(numel(ND)<=numel(MDEL1))
                ND=1:n;
                MDEL1=MDEL;
            end
            ND(MDEL1)=[];
            nd=ND(randi(length(ND),1,1));
            ud=(rand<0.5);
            fprintf('Randomly decimated at count_n=%d\n',count_n);
        end
        UD(nd)=ud;
        MDEL=[MDEL nd];
        parities_error_sequence=ParityError(UD,MDEL,Vcon,Ccon1);
        if(sum(parities_error_sequence)-prev_parities_error>0)
            if(RAND_DECIM)
                fprintf('Randomly decimated in this iter\n');
            end
            fprintf('count_n=%d,parities in error=%d\n',count_n,sum(parities_error_sequence));
        end
        prev_parities_error=sum(parities_error_sequence);
        count_n=count_n-1;
        yr=Updatey(yr,Vcon,Ccon2,nd,ud);
        if(mod(count_n,100)==0)
            fprintf('Finished %d/%d code-bits, par in err/total parities=%d/%d\n',n-count_n,n,sum(parities_error_sequence),m1);
        end
    end
    
    E_check_code(Vcon)=repmat(UD,1,dv);
    parities_in_error=sum(mod(sum(E_check_code(Ccon1),2),2));
    if parities_in_error==0
        fprintf('Encoding successful\n');
    else
        fprintf('Not Encoded Correctly!, parities in error=%d\n',parities_in_error);
    end
    eps_ch=parities_in_error*dc1/n;
    Distortion(j)=sum(yr)/m;
    DistProfile(:,j)=sum(reshape(yr,[],L),1)'/ms;

    fprintf('Distortion for the curr. seq is %f, avg is %f \n',Distortion(j), mean(Distortion(1:j)) );
%     fprintf('Shannon Distortion is %f\n',dist_sha);
 end
