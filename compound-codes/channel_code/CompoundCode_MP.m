clear; 
%clc;

% For the LDPC Part
dv1=1; % Variable-node degree
dc1=4; % Check-node degree

%%%%%%%%%%
% This requires n=2*m1
%%%%%%%%%%

% For the LDGM Part
dv2=9; % Variable-node degree
dc2=5; % Check-node degree

%%%%%%%%%%%
% LDGM Part: L=x^(dv2); R=tx+(1-t)x^(dc2);
% n * dv2 = m * (t + (1-t)dc2);
% This requires m*t to be an integer
%%%%%%%%%%%
t=0.3;
% Choosing t=0.3; This gives 45*n=19*m with dv2=9 and dc2=5;

c=1e4;
m=45*(c); % Block Length / # of check-nodes in LDGM
n=round(m*(t+(1-t)*dc2)/dv2); % # of variable-nodes
m1=n*dv1/dc1; % # of check-nodes in LDPC
iter=200; % Max iterations of Message Passing algorithm
M=10; % # of codewords for averaging

e1=dv1*n; % # of edges in LDPC part
e2=dv2*n; % # of edges in the LDGM part
e=e1+e2; % Total # of edges

% Extra edge to handle irregularity in LDGM part
E=zeros(1,e+1); % This is the main object in the sim

Vcon=[reshape(1:e1,dv1,[])' reshape(e1+1:e,dv2,[])']; % variable-node connections
E_perm=[randperm(e1) e1+randperm(e2)];
Ccon1=reshape(E_perm(1:e1),[],dc1); % check-node connections in LDPC

% Handle irregularity at the check-nodes in LDGM
Ccon2=zeros(m,dc2); % check-node connections in LDGM
Ccon2(1:(t*m),1)=E_perm((e1+1):(e1+t*m)); Ccon2(1:(t*m),2:end)=e+1;
Ccon2((t*m+1):end,:)=reshape(E_perm((e1+t*m+1):end),[],dc2);

E_code=zeros(1,e+1); % Object to determine the codeword given code-bits

% errV=0.01:0.001:0.13;
errV=0.09:0.003:0.12;
Pb=zeros(size(errV)); % Bit error probability

for k=1:1%length(errV)
    fprintf('Simulating for error prob of %f\n',errV(k));
    err=errV(k);
    llr0=log((1-err)/err);
    Eold=zeros(size(E));
    Pb_M=zeros(numel(errV),M);
    for j=1:M
        tic;
    y=ones(m,1); y(rand(m,1)<err)=-1;
        for i=0:iter
            if(i==0)
                % Routine for first iteration
                E(Ccon1)=0;
                E(Ccon2(1:(m*t),:))=repmat(llr0*y(1:(m*t)),1,dc2);
                E(Ccon2((m*t+1):end,:))=0;
                continue;
            end
            %Variable-node update
            E(Vcon)=repmat(sum(E(Vcon),2),1,dv1+dv2)-E(Vcon);
            E(E>20)=20; E(E<-20)=-20;
            E(end)=0;   
            
            % Checking the error assuming hard decision decoding
            V=sum(E(Vcon),2)<0;
            E_code(Vcon)=repmat(V,1,dv1+dv2); E_code(end)=0;
            checks_in_error=sum(mod(sum(E_code(Ccon1),2),2));
            bits_in_error=sum(mod(sum(E_code(Ccon2),2),2));
            
            fprintf('After %d iter, (bits,checks) in error: (%d,%d)\n',i,bits_in_error,checks_in_error);
            if(max(abs(Eold-E))<1e-9)
                fprintf('Messages Converged! Stopping...\n');
                break;
            end
            Eold=E;
            
            % Check-node update
            E(end)=Inf; E(E==0)=eps;
            E(Ccon1)=repmat(prod(tanh(E(Ccon1)/2),2),1,dc1)./tanh(E(Ccon1)/2);
            E(Ccon2)=repmat(tanh(llr0/2)*prod(tanh(E(Ccon2)/2),2).*y,1,dc2)./tanh(E(Ccon2)/2);
            E=2*atanh(E);
            E(E>20)=20; E(E<-20)=-20;
        end
        fprintf('j/M=%d/%d, Pb_M=%f\n',j,M,bits_in_error/m);
        Pb_M(k,j)=bits_in_error;
    toc
    end
    Pb(k)=mean(Pb_M(k,:));
end