function xhat=BP_2GMAC(y,channel,Vcon,Ccon,eMax,maxIters,Parity,shiftPattern,sig)
% dv and dc stand for Max Variable and Check node degrees respectively.
% This is intended for a Irregular degree profile too.
% iter=Max iterations of mess. pass.
global llr_max
[N,dv]=size(Vcon);
m=length(Ccon(:,1));
llr_max=20;
Ccon(Ccon>eMax)=eMax+1;
invshiftPattern(shiftPattern)=1:N;
if 2*m== numel(Parity)
    if m~= length(Parity(:,1))
        Parity=Parity';
    end
    U=numel(Parity(1,:));
else
    error('Length of Parity of coset does not match with Ccon');
end

E1=zeros(eMax,1); E1(eMax+1)=Inf; E2=E1;
V1_btoc=zeros(N,1);      V2_btoc=V1_btoc;
opLLR=zeros(length(V1_btoc),2);

%-----Message Passing------------------------------------------
for i=0:maxIters
% Check-node update            
    E1=2*atanh(extr_prod(E1,Ccon,20*sign(1-2*Parity(:,1))));           
    E2=2*atanh(extr_prod(E2,Ccon,20*sign(1-2*Parity(:,2))));           
    
    E1(eMax+1)=0;  E2(eMax+1)=0;    
    
    %Update at MAC node (real sum)
   llr_out=LLR_real_adder(reshape(y,[],1),[V1_btoc V2_btoc(shiftPattern)],channel,U,sig) ;
   if ~isempty(find(isnan(llr_out),1))
        break
    end
   V1_ctob=llr_out(:,1);  V2_ctob=llr_out(invshiftPattern,2);
    
% Variable-node update
    V1_fromChecks=sum(E1(Vcon),2);   V2_fromChecks=sum(E2(Vcon),2);
    E1(Vcon)=repmat(V1_fromChecks+V1_ctob,1,dv)-E1(Vcon);
    E2(Vcon)=repmat(V2_fromChecks+V2_ctob,1,dv)-E2(Vcon);
    
    E1(E1>20)=20; E1(E1<-20)=-20; E1(abs(E1)<1e-5)=0;
    E2(E2>20)=20; E2(E2<-20)=-20; E2(abs(E2)<1e-5)=0;
    E1(eMax+1)=Inf;  E2(eMax+1)=Inf;
    
    %Update at bit node (message to MAC node)
    V1_btoc=V1_fromChecks;  V2_btoc=V2_fromChecks;
    
    if mod(i,10)==0
        opLLR_prev=opLLR;
        opLLR=[V1_fromChecks+V1_ctob V2_fromChecks+V2_ctob]; 
    end
    if max(max(abs(opLLR_prev-opLLR)))<0.001
       break 
    end   
    if ~isempty(find(isnan(E1),1)) || ~isempty(find(isnan(E2),1))
        break
    end
end
xhat=sign([V1_fromChecks+V1_ctob V2_fromChecks+V2_ctob])';
end