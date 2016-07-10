clear; 
% clc;

% For the LDPC Part
dv1=3; % Variable-node degree
dc1=6; % Check-node degree

%%%%%%%%%%%%%%%%%%%%%%%
% This requires ns*dv1=m1s*dc1
%%%%%%%%%%%%%%%%%%%%%%%

% For the LDGM Part
dv2=3; % Variable-node degree
dc2=3; % Check-node degree

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LDGM Part: L=x^(dv2); R=tx+(1-t)x^(dc2);
% n * dv2 = m * (t + (1-t)dc2);
% This requires m*t to be an integer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=0.2;
% Choosing t=0.3; This gives 45*n=19*m with dv2=9 and dc2=5;

c=48;
ms=300*(c); % # of check-nodes in LDGM in Single System
ns=round(((t+(1-t)*dc2)/dv2)*ms); % # of variable-nodes in a Single System
% ns=149*(c); % # of variable-nodes in a Single System
m1s=ns*dv1/dc1; % # of check-nodes in LDPC in a Single System

rate=(1-dv1/dc1)*(ns/ms);
eps_sha=finverse(@(x)(BinaryEntropy(x)),1-rate,0,0.5);
fprintf('rate=%f,e_shannon is %f\n',rate,eps_sha);

L=15; % Chain length of the SC system
w=3; % coupling parameter

n=ns*L; % Total # of variable-nodes
m=ms*L; % Total # of check-nodes in LDGM / Block length
m1=m1s*L; % Total # of check-nodes in LDPC

e1s=dv1*ns; % # of edges in LDPC part in single system
e2s=dv2*ns; % # of edges in the LDGM part in single system
es=e1s+e2s; % Total # of edges in single system

e1=e1s*L;
e2=e2s*L;
e=es*L; % Total # of edges

% Extra edge to handle irregularity in LDGM part
E=zeros(1,e+1); % This is the main object in the sim

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
    Ccon2(((i-1)*ms+1) : ((i-1)*ms+ms*t),1)=Etmp2(1:ms*t);
    Ccon2(((i-1)*ms+1) : ((i-1)*ms+ms*t),2:end)=e+1;
    Ccon2(((i-1)*ms+ms*t+1) : (i*ms),:)=reshape(Etmp2((ms*t+1):end),[],dc2);
end

E_code=zeros(1,e+1); % Object to determine the codeword given code-bits

iter=200; % Max iterations of Message Passing algorithm
M=10; % # of codewords for averaging

% errV=0.15:0.01:0.19;
errV=0.128:0.001:0.205;
% errV=0.08;
Pb=zeros(size(errV)); % Bit error probability
Pb_M=zeros(numel(errV),M);

for k=1:length(errV)
fprintf('Simulating for error prob of %f\n',errV(k));
err=errV(k);
llr0=log((1-err)/err);
for j=1:M
E=zeros(1,e+1);
Eold=zeros(size(E));
y=ones(m,1); y(rand(m,1)<err)=-1;
tic; 
for i=1:iter
%     if(i==0)
%     % Routine for first iteration
%     E(Ccon1)=0;
%     for z=1:L
%         E(Ccon2(((z-1)*ms+1):((z-1)*ms+ms*t),:))=repmat(llr0*y(((z-1)*ms+1):((z-1)*ms+ms*t)),1,dc2);
%         E(Ccon2(((z-1)*ms+ms*t+1) : (z*ms),:))=0;
%     end
%     continue;
%     end
    %Variable-node update
    E(Vcon)=repmat(sum(E(Vcon),2),1,dv1+dv2)-E(Vcon);
    E(Vcon(1:((w-1)*ns),:))=20; % boundary condition in a SC system
    E(E>20)=20; E(E<-20)=-20;
    E(end)=0;

    % Check-node update
    E(end)=Inf;
    E(Ccon1)=GenProd(E(Ccon1),1/2);
    E(end)=Inf;
    E(Ccon2)=repmat(tanh(llr0/2).*y,1,dc2).*GenProd(E(Ccon2),1/2);
    E=2*atanh(E);
    E(E>20)=20; E(E<-20)=-20;
    E(end)=0;

    % Checking the error assuming hard decision decoding 
    V=sum(E(Vcon),2)<0;
    E_code(Vcon)=repmat(V,1,dv1+dv2); E_code(end)=0;
    checks_in_error=sum(mod(sum(E_code(Ccon1),2),2));
    bits_in_error=sum(mod(sum(E_code(Ccon2),2),2));
    if(mod(i,5)==0)
        fprintf('After %d iter, (bits,checks) in error: (%d,%d)\n',i,bits_in_error,checks_in_error);
    end
    if(max(abs(Eold(1:e)-E(1:e)))<1e-9)
        fprintf('Messages Converged! Stopping...\n');
        break;
    end
    Eold=E;
end
toc;
fprintf('Error=%f\n',bits_in_error/m);
Pb_M(k,j)=bits_in_error/m;
end
Pb(k)=mean(Pb_M(k,:));
end
hold on;
plot(errV,Pb);
