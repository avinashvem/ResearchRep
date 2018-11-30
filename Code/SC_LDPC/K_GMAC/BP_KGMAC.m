function xhat=BP_KGMAC(y,Vcon,Ccon,eMax,maxIters,Parity,sig)
global interl deinterl
U=50;


%{
N-#bit nodes in the grpah, m-#check nodes
K-#of users in GMAC.
y=\sum_{i=1,..,K} pi_i(x_i) +z where pi_i is the ith interleaver-interl(:,i)

% y-Nx1
% eMax-Max. # edges in the graph
% Parity-mxU 
% interl- NxK where m #check nodes, K-#of users in GMAC
%}
[n,dv]=size(Vcon); [m,~]=size(Ccon);
N=numel(y);
Ccon(Ccon>eMax)=eMax+1;
K=numel(Parity(1,:));

G=exp(-((repmat(y,1,2*U+1)- repmat(-U:U,N,1)).^2)/2/sig^2);
if m*K== numel(Parity)
    if m~= length(Parity(:,1))
        Parity=Parity';
    end
else
    error('Length of Parity of coset does not match with Ccon');
end

E1=zeros(eMax,1); E1(eMax+1)=Inf; E2=E1;
E=repmat(E1,1,K);

E_btog=[zeros(n,K);1i*ones(N-n,K)];  opLLR=E_btog(1:n,:);


%--------Message Passing------------------------------------
for iter=0:maxIters
% Check-node update
    for u=1:K
        E(:,u)=2*atanh(extr_prod(E(:,u),Ccon,20*sign(1-2*Parity(:,u))));           
    end
    %{
    E1=2*atanh(extr_prod(E1,Ccon,20*sign(1-2*Parity(:,1))));           
    E2=2*atanh(extr_prod(E2,Ccon,20*sign(1-2*Parity(:,2))));           
    %}
    E(eMax+1,:)=0;
% Update at MAC node (real sum) 
   for u=1:K
      perm_E_btog(:,u)=E_btog(interl(:,u),u);
   end
   llr_out=LLR_real_adder(perm_E_btog,y,sig,G);
   
   if ~isempty(find(isnan(llr_out),1))
        break
   end
   
   for u=1:K
      E_gtob(:,u)=llr_out(deinterl(1:n,u),u);
   end
    
% Variable-node update (message to check node)
    for u=1:K
        Eu=E(:,u);  E_MAC(:,u)=sum(Eu(Vcon),2);  
        
        Eu(Vcon)=repmat(E_MAC(:,u)+E_gtob(:,u),1,dv)-Eu(Vcon);
    
        Eu(Eu>20)=20; Eu(Eu<-20)=-20; Eu(abs(Eu)<1e-5)=0;
        Eu(eMax+1)=Inf;  
        E(:,u)=Eu;
        
% Update at bit node (message to MAC node)
    end
    E_btog=[E_MAC;1i*ones(N-n,K)];
    
    if mod(iter,5)==0
        opLLR_prev=opLLR;
        opLLR=E_MAC+E_gtob; 
    end
    if max(max(abs(opLLR_prev-opLLR)))<0.001
       break 
    end   
    if ~isempty(find(isnan(E1),1)) || ~isempty(find(isnan(E2),1))
        break
    end
end
xhat=[sign(real(E_MAC+E_gtob));zeros(N-n,K)];

for u=1:K
  xhat(:,u)=xhat(interl(:,u),u);
end
 end