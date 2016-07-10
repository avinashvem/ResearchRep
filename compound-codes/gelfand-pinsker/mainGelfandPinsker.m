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
% This gives 5*n=2*m with dv2=6 and dc2=3;

% Random qth fraction of LDPC check nodes 
% in each system are not part of channel code.
q=0.5;

c=12;
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

isFix4Cycles=0;
while(~isFix4Cycles)
    isFix2Cycles=0;
    while(~isFix2Cycles)
        [Vcon Ccon1 Ccon2 isFix2Cycles]=fix2cycles(Vcon,Ccon1,Ccon2,dv1,dv2);
    end
    [Vcon Ccon1 Ccon2 isFix4Cycles]=fix4cycles(Vcon,Ccon1,Ccon2,dv1,dv2);
end

Vcon=Vcon(1:(L*ns),:);
E_code=zeros(1,e+1);
E_code(Vcon)=ones(size(Vcon));
Ccon2(E_code(Ccon2)==0)=e+2;
n=ns*L;

% I denotes which LDPC check-nodes are active in channel code
I=ones(m1,1); I(rand(m1,1)<q)=0; I(1:((w-1)*m1s))=1; I(((ldpc_len-w+1)*m1s+1):(ldpc_len*m1s))=1; I=logical(I);
Ccon1_sub=Ccon1(I,:); MSG_I=~I;
msg_size=sum(MSG_I);

%%%% Delete unused Vcon sockets
E_code=zeros(1,e+2);
E_code(Ccon1)=ones(size(Ccon1)); E_code(Ccon2)=ones(size(Ccon2));
Vcon(E_code(Vcon)==0)=e+1;

%%%% Vcon_sub denotes Vcon connections for subcode
Vcon_sub=Vcon;
E_code=zeros(1,e+2);
E_code(Ccon1_sub)=ones(size(Ccon1_sub)); E_code(Ccon2)=ones(size(Ccon2));
Vcon_sub(E_code(Vcon_sub)==0)=e+1;

M=3; % # of source sequences to average the distortion

rate_chan=(ns-(1-q)*m1s)/ms; 
chan_eps_sha=finverse(@(x)(BinaryEntropy(x)),1-rate_chan,0,0.5);
fprintf('Channel coding rate=%f, Shannon Th is %f\n',rate_chan,chan_eps_sha);

rate_dist=(ns-m1s)/ms; %Total Rate
dist_sha=finverse(@(x)(BinaryEntropy(x)),1-rate_dist,0,0.5);
fprintf('Total rate=%f, Shannon Distortion is %f\n',rate_dist,dist_sha);

Dist=zeros(M,1);
DistProfile=zeros(M,L);
UD_Store=zeros(M,n);
YR_Store=zeros(M,m);
S_Store=zeros(M,m);
Msg_Store=zeros(M,msg_size);

ERR_out=zeros(M,1);

for i=1:M
    % Random side-information
    S=zeros(m,1); S(rand(m,1)<0.5)=1;
    S(1:((w-1)*ms))=0;
    Msg=zeros(msg_size,1); Msg(rand(msg_size,1)<0.5)=1;
    BinIndexOrg=zeros(m1,1); BinIndexOrg(MSG_I)=Msg;
    [Dist_tmp DistProfile_tmp UD yr]=CompressSource(Vcon,Ccon1,Ccon2,S,MSG_I,Msg);
    Dist(i)=Dist_tmp;
    DistProfile(i,:)=DistProfile_tmp;
    UD_Store(i,:)=UD;
    YR_Store(i,:)=yr;
    S_Store(i,:)=S;
    Msg_Store(i,:)=Msg;    
    disp(Dist_tmp); disp(DistProfile_tmp);
    save 3_6_10_5_t=0.01_beta=0.7.mat
end

return;

err=0.9*chan_eps_sha;
llr0=log((1-err)/err);
E_code=zeros(1,e+1);
E_code(Vcon)=repmat(UD,1,dv1+dv2);
S_hat=mod(sum(E_code(Ccon2),2),2);

W=zeros(m,1); W(rand(m,1)<err)=1;
Y=xor(S_hat,W);
UD_Dec=MP_CompoundSC(Vcon_sub,Ccon1_sub,Ccon2,Y,llr0);

E_code=zeros(1,e+1);
E_code(Vcon)=repmat(UD_Dec,1,dv1+dv2);
BinIndex=mod(sum(E_code(Ccon1),2),2);
Msg_out=BinIndex(MSG_I);
ERR_out(i)=mean(xor(Msg,Msg_out));
disp(ERR_out(i));
