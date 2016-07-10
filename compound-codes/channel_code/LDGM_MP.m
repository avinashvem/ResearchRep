clear; clc;

dv=10; % Variable-node degree
dc=5; % Check-node degree

%%%%%%%%%%%
% LDGM Part: L=x^(dv); R=tx+(1-t)x^(dc);
% n * dv = m * (t + (1-t)dc);
% This requires m*t to be an integer
%%%%%%%%%%%
t=0.3;
% Choosing t=0.3; This gives 45*n=19*m with dv=9 and dc=5;

m=50*(2e4); % Block Length / # of check-nodes in LDGM
n=19*(2e4); % Number of variable-nodes
iter=200; % Max iterations of Message Passing algorithm
M=1; %Number of codewords for averaging

e=dv*n; % Number of edges in the LDGM part
% Extra edge to handle irregularity
E=zeros(1,e+1); % This is the main object in the simulation
E_code=zeros(size(E)); % Given the code-bits, this finds the codeword
Vcon=reshape(1:e,[],dv);

% Handle irregularity at the check-nodes
E_perm=randperm(e);
Ccon=zeros(m,dc);
Ccon(1:(t*m),1)=E_perm(1:(t*m)); Ccon(1:(t*m),2:end)=e+1;
Ccon((t*m+1):end,:)=reshape(E_perm((t*m+1):end),[],dc);

errV=0.1:0.001:0.115;
% errV=0.090;
Pb=zeros(size(errV)); % Bit error probability

for k=1:length(errV)
    fprintf('Simulating for error prob of %f\n',errV(k));
    err=errV(k);
    llr0=log((1-err)/err);
    Eold=zeros(size(E));
    Pb_M=zeros(numel(errV),M);
    for j=1:M
    y=ones(m,1); y(rand(m,1)<err)=-1;
        for i=0:iter
            if(i==0)
                % Routine for first iteration
                E(Ccon(1:(m*t),:))=repmat(llr0*y(1:(m*t)),1,dc);
                E(Ccon((m*t+1):end,:))=0;
                continue;
            end
            %Variable-node update
            E(Vcon)=repmat(sum(E(Vcon),2),1,dv)-E(Vcon);
            E(E>20)=20; E(E<-20)=-20;
            
            % Checking the error assuming hard decision decoding 
            V=sum(E(Vcon),2)<0;
            E_code(Vcon)=repmat(V,1,dv); E_code(end)=0;
            bits_in_error=sum(mod(sum(E_code(Ccon),2),2));
            
            fprintf('After %d iter, bits in err: %d\n',i,bits_in_error);    
            if(max(abs(Eold-E))<1e-7)
                fprintf('Messages Converged! Stopping...\n');
                break;
            end
            Eold=E;
            % Check-node update
            E(e+1)=Inf; E(E==0)=eps;
            E(Ccon)=repmat(tanh(llr0/2)*prod(tanh(E(Ccon)/2),2).*y,1,dc)./tanh(E(Ccon)/2);
            E(Ccon)=2*atanh(E(Ccon));
            E(E>20)=20; E(E<-20)=-20;
        end
        Pb_M(k,j)=bits_in_error/m;
    end
    Pb(k)=mean(Pb_M(k,:));
end
hold on;
plot(errV,Pb,'-r');
