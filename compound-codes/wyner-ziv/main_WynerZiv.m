clear; clc;
global L w;
% For the LDPC part
dv1=3; % Variable-node degree
dc1=6; % Check-node degree

% For the LDGM Part
dv2=6; % Variable-node degree
dc2=3; % Check-node degree
dv=dv1+dv2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LDGM Part: L=x^(dv2); R=tx+(1-t)x^(dc2);
% n * dv2 = m * (t + (1-t)dc2);
% This requires m*t to be an integer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=0.01;
% Choosing t=0.3; This gives 45*n=19*m with dv2=9 and dc2=5;

c=4;
ms=300*(c); % # of check-nodes in LDGM in Single System
ns=149*(c); % # of variable-nodes in a Single System
m1s=ns*dv1/dc1; % # of check-nodes in LDPC in a Single System
e1s=ns*dv1; e2s=ns*dv2; es=ns*dv;

L=20; % Chain length of the SC system
w=4; % coupling parameter
Lw=L+w-1; % code-bits system length
ldpc_start=w; ldpc_end=L; 
ldpc_len=ldpc_end-ldpc_start+1; % ldpc-check system length
m=ms*L;
n=ns*Lw;
m1=m1s*ldpc_len;

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
    
    Ccon2(((i-1)*ms+1) : ((i-1)*ms+ms*t),1)=Etmp2(1:ms*t);
    Ccon2(((i-1)*ms+1) : ((i-1)*ms+ms*t),2:end)=e+1;
    Ccon2(((i-1)*ms+ms*t+1) : (i*ms),:)=reshape(Etmp2((ms*t+1):end),[],dc2);
end

%%%% Delete unused Vcon sockets
E_code=zeros(1,e+1);
E_code(Ccon1)=ones(size(Ccon1)); E_code(Ccon2)=ones(size(Ccon2));
Vcon(E_code(Vcon)==0)=e+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M=10; % # of source sequences to average the distortion

rate_ldgm=ns/ms; 
dist_sha=finverse(@(x)(BinaryEntropy(x)),1-rate_ldgm,0,0.5);
fprintf('ldgm_rate=%f, Shannon distortion is %f\n',rate_ldgm,dist_sha);

rate=(ns-m1s)/ms;
chan_eps_sha=finverse(@(x)(BinaryEntropy(x)),1-rate,0,0.5);
fprintf('rate=%f, Channel Shannon Th is %f\n',rate,chan_eps_sha);

Dist=zeros(M,1);
DistOut=zeros(M,1);
DistProfile=zeros(M,L);
DistProfile_final=zeros(M,L);

UD_Store=zeros(M,n);
y_Store=zeros(M,m);
yr_Store=zeros(M,m);

for i=1:M
    % Bernoulli(1/2) symmetric source
    y=zeros(m,1); y(rand(m,1)<0.5)=1; % current instance of source seq.
    y(1:((w-1)*ms))=0;
    [Dist_tmp DistProfile_tmp UD yr]=CompressSource(Vcon(:,(dv1+1):end),Ccon2,y,ns);
    Dist(i)=Dist_tmp; DistProfile(i,:)=DistProfile_tmp;
    save 3_6_6_3_t=0.01.mat
    disp(Dist_tmp); disp(DistProfile_tmp);
    UD_Store(i,:)=UD; y_Store(i,:)=y; yr_Store(i,:)=yr;
end

return;

% Post-Processing;

err=0.92*chan_eps_sha;
err_eff=(err-Dist(i))/(1-2*Dist(i));

E_code=zeros(1,e+1);
E_code(Vcon)=repmat(UD,1,dv1+dv2);
BinIndex=mod(sum(E_code(Ccon1),2),2);

%     y_test=mod(sum(E_code(Ccon2),2),2);

y_side_inf=zeros(size(y)); y_side_inf(rand(m,1)<err_eff)=1;
y_side_inf=xor(y,y_side_inf); y_side_inf(1:((w-1)*ms))=0;

UD_Dec=ChanDecBinIndex(Vcon,Ccon1,Ccon2,BinIndex,err,y_side_inf);
E_code=zeros(1,e+1);
E_code(Vcon)=repmat(UD_Dec,1,dv1+dv2);

y_decoded=mod(sum(E_code(Ccon2),2),2); y_decoded(1:((w-1)*ms))=0;
DistProfile_final(i,:)=sum(reshape(xor(y_decoded,y),[],L),1)/(m/L);
disp(DistProfile_final(i,:));
DistOut(i)=mean(DistProfile_final(i,w+2:L-w));
disp(DistOut(i));
