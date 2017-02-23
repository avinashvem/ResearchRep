function xhat=BP_two_user_MAC_2(y,channel,Vcon,Ccon,eMax,maxIters,Parity,sig)
% This is exactly identical to the BP_two_user_MAC except for the graph
% relizatons for the two users. We use two independent Tanner graphs from
% the same ensemble for the two users.

% dv,dc are Max Variable and Check node degrees respectively. Works for Irregular d.ds.
% y=x1+x2(shiftPattern)+n(0,sig^2) hwere shiftPattern is a random shift at one of the users.

global llr_max
N=length(Vcon(:,1))/2; dv=length(Vcon(1,:));
m=length(Ccon(:,1))/2;
Ccon(Ccon>eMax)=eMax+1;

Vcon1=Vcon(1:N,:); Vcon2=Vcon(N+(1:N),:);
Ccon1=Ccon(1:m,:); Ccon2=Ccon(m+(1:m),:);

llr_max=20;
if 2*m== numel(Parity)
    if m~= length(Parity(:,1))
        Parity=Parity';
    end
    U=numel(Parity(1,:));
else
    error('Length of Parity of coset does not match with Ccon');
end
%{
if U==2
   if strcmp(channel,'noiseless')
       llr1=zeros(numel(y),1);        
       pos= y==2; neg= y==-2;
       llr1(pos)=llr_max; llr1(neg)=-llr_max; 
       V1Chan=llr1; V2Chan=V1Chan(invshiftPattern);
       nonZeroIdx1= V1Chan~=0;     nonZeroIdx2= V2Chan~=0;
   end
end
%}

E1=zeros(eMax,1); E1(eMax+1)=Inf; E2=E1;
V1_btoc=zeros(N,1);      V2_btoc=V1_btoc;
%% Message Passing
for i=1:maxIters
    E1_prev=E1; E2_prev=E2;
% Check-node update            
    E1=2*atanh(extr_prod(E1,Ccon1,20*sign(1-2*Parity(:,1))));           
    E2=2*atanh(extr_prod(E2,Ccon2,20*sign(1-2*Parity(:,2))));           
    
    E1(E1>20)=20; E1(E1<-20)=-20;
    E2(E2>20)=20; E2(E2<-20)=-20;
    E1(eMax+1)=0;  E2(eMax+1)=0;    
    
% Update at MAC node (real sum)
    llr_out=LLR_real_adder(reshape(y,[],1),[V1_btoc V2_btoc],channel,U,sig) ;
    V1_ctob=llr_out(:,1);  V2_ctob=llr_out(:,2);
    
% Variable-node update
    V1_fromChecks=sum(E1(Vcon1),2);   V2_fromChecks=sum(E2(Vcon2),2);
    E1(Vcon1)=repmat(V1_fromChecks+V1_ctob,1,dv)-E1(Vcon1);
    E2(Vcon2)=repmat(V2_fromChecks+V2_ctob,1,dv)-E2(Vcon2);
    
    E1(E1>20)=20; E1(E1<-20)=-20; E1(abs(E1)<1e-5)=0;
    E2(E2>20)=20; E2(E2<-20)=-20; E2(abs(E2)<1e-5)=0;
    E1(eMax+1)=Inf;  E2(eMax+1)=Inf;
    
    %Update at bit node (message to MAC node)
    V1_btoc=V1_fromChecks;  V2_btoc=V2_fromChecks;
    
    if max(E1_prev-E1)<0.00001 && max(E2_prev-E2)<0.00001
       break 
    end   
end
xhat=sign([V1_fromChecks+V1_ctob V2_fromChecks+V2_ctob])';
end
