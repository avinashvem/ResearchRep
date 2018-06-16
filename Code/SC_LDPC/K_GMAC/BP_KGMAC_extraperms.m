function xhat=BP_KGMAC_extraperms(y,Vcon,Ccon,eMax,maxIters,Parity,sig,interl,deinterl)
U=50;
llr_max=20;
llr_min=1e-7;

%{
N-#bit nodes in the grpah, m-#check nodes
Ka-#of users in GMAC.
Kb-#of ineterleavers given as input
y=\sum_{i=1,..,Ka} pi_i(x_i) +z where pi_i is the ith interleaver-interl(:,i)

% y-Nx1
% eMax-Max. # edges in the graph
% Parity-mx(Ka+Kb) where m #check nodes
% interl- Nx(Ka+Kb) 
%}
[n,dv]=size(Vcon); [m,~]=size(Ccon);
N=numel(y);
Ccon(Ccon>eMax)=eMax+1;
Ka=numel(Parity(1,:));
Kb=numel(interl(1,:))-Ka;
Parity=[Parity zeros(m,Kb)];
K=Ka+Kb;

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
        E(:,u)=2*atanh(extr_prod(E(:,u),Ccon,llr_max*sign(1-2*Parity(:,u))));           
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
    
        Eu(Eu>llr_max)=llr_max; Eu(Eu<-llr_max)=-llr_max; Eu(abs(Eu)<llr_min)=0;
        Eu(eMax+1)=Inf;  
        E(:,u)=Eu;
        
% Update at bit node (message to MAC node)
    end
    E_btog=[E_MAC;1i*ones(N-n,K)];
%    E_MAC(E_MAC>llr_max)=llr_max; E_MAC(E_MAC<-llr_max)=-llr_max;
%   E_gtob(E_gtob>llr_max)=llr_max; E_gtob(E_gtob<-llr_max)=-llr_max;
    

    if mod(iter,5)==0
        opLLR_prev=opLLR;
        opLLR=E_MAC+E_gtob; 
        opLLR(opLLR>llr_max)=llr_max; opLLR(opLLR<-llr_max)=-llr_max;
        max_change=max(max(abs(opLLR_prev-opLLR)));
        fprintf('%3.2f,',max_change);
     if max_change<0.001
       break 
     end   
    end
    
    if ~isempty(find(isnan(E1),1)) || ~isempty(find(isnan(E2),1))
        break
    end
end

abs_llr=sum(abs(real(E_MAC+E_gtob)))/n;
idx=zeros(1,Kb);
for i=1:Kb
    [~,idx(i)]=min(abs_llr);
    abs_llr(idx(i))=50;
end
rem_idx=setdiff(1:K,idx);

xhat=[sign(real(E_MAC+E_gtob));zeros(N-n,K)];
for i=1:Ka
    u=rem_idx(i);
  xhat(:,u)=xhat(interl(:,u),u);
end
 end