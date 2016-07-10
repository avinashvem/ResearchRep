clear; 
% clc;

% For the LDPC Part
dv1=3; % Variable-node degree
dc1=6; % Check-node degree

% For the LDGM Part
dv2=6; % Variable-node degree
dc2=3; % Check-node degree
dv=dv1+dv2;
%%%%%%%%%%%

rate=(dc2/dv2)*(1-(dv1/dc1)); 
dist_sha=finverse(@(x)(BinaryEntropy(x)),1-rate,0,0.5);
fprintf('rate=%f,Shannon distortion is %f\n',rate,dist_sha);

c=4e3;
m=dv2*(c); % Block Length / # of check-nodes in LDGM
n=dc2*(c); % # of variable-nodes
m1=n*(dv1/dc1); % # of check-nodes in LDPC

e1=dv1*n; % # of edges in LDPC part
e2=dv2*n; % # of edges in LDGM part
e=e1+e2; % Total # of edges

% Extra edge to handle irregularity in LDGM part
E=zeros(1,e); % This is the main object in the sim
E_check_code=zeros(1,e); % Object to check if the code-bits satisfy LDPC checks

Vcon=[reshape(1:e1,dv1,[])' reshape(e1+1:e,dv2,[])']; % variable-node connections
Ccon1=reshape(randperm(e1),[],dc1); % check-node connections in LDPC

Ccon2=reshape(e1+randperm(e2),[],dc2);

M=1; % # of source sequences to average the distortion

beta=1; % inverse temperature
T=10; % Number of BP iterations before a forced decimation
alp=4.25; % threshold for decimation

% Messages on the edges

Distortion=zeros(1,M);
tic;
for j=1:M
    UD=zeros(n,1);
    % Bernoulli symmetric source
    y=zeros(m,1); y(rand(m,1)<0.5)=1; % current instance of source seq.
    yr=y;
    E=zeros(1,e);
    count_n=n;
    RAND_DECIM=-1*ones(2,n);
    BIAS_DECIM=zeros(n,dv);
    MDEL=[]; % Variable bits deleted
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
            I=find(abs(bt)==B);
            nd=I(randi(length(I),1,1));
            ud=bt(nd)<0;
            RAND_DECIM(1,nd)=0; RAND_DECIM(2,nd)=B;
            BIAS_DECIM(nd,:)=E(Vcon(nd,:));
            assert(all(MDEL-nd));
        else
            ND=1:n; ND(MDEL)=[];
            nd=ND(randi(length(ND),1,1));
            RAND_DECIM(1,nd)=1; RAND_DECIM(2,nd)=0;
            BIAS_DECIM(nd,:)=E(Vcon(nd,:));
            ud=(rand<0.5);
        end
        UD(nd)=ud;
        MDEL=[MDEL nd];
        parities_error_sequence=ParityError(UD,MDEL,Vcon,Ccon1);
        if(sum(parities_error_sequence))
%             fprintf('count_n=%d,parities in error=%d\n',count_n,sum(parities_error_sequence));
        end
        count_n=count_n-1;
        yr=Updatey(yr,Vcon,Ccon2,nd,ud);
        if(mod(count_n,100)==0)
            fprintf('Finished %d/%d code-bits, par in err=%d\n',n-count_n,n,sum(parities_error_sequence));
        end
    end
    
    E_check_code(Vcon)=repmat(UD,1,dv);
    parities_in_error=sum(mod(sum(E_check_code(Ccon1),2),2));
    if parities_in_error==0
        fprintf('Encoding successful\n');
    else
        fprintf('Not Encoded Correctly!, parities in error=%d\n',parities_in_error);
    end
    Distortion(j)=sum(yr)/m;
    fprintf('Distortion for the curr. seq is %f, avg is %f \n',Distortion(j), mean(Distortion(1:j)) );
    fprintf('Shannon Distortion is %f\n',dist_sha);
end
toc;