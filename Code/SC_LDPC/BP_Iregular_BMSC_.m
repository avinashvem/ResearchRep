function op=BP_Iregular_BMSC(V,Vcon,Ccon,eMax,Parity)
% This is intended for a Irregular degree profiles.
% Assumes all zero codeword especially break condition for Message passing to stop.
% iter signifies Max iterations of Message Passing algorithm
e=eMax;
[N,dv]=size(Vcon);
[m, dc]=size(Ccon);
if m== numel(Parity)
    if m~= length(Parity(:,1))
        Parity=Parity';
    end
else
    error('Parity of coset does not match with Ccon');
end
if numel(V(1,:))~=1
   V=V'; 
end

E=zeros(e,1);
E(e+1)=Inf;
Ccon(Ccon>e)=e+1;
E_prev=zeros(length(E),1);
E(Vcon)=repmat(V,1,dv);
%------Message Passing----------------------

 % Check-node update            
    E=extr_prod(E,Ccon,20*sign(1-2*Parity));           
    E=2*atanh(E);
    E(E>20)=20; E(E<-20)=-20;
    E(e+1)=0;
  % Variable-node update
%     E(Vcon)=repmat(sum(E(Vcon),2)+V,1,dv)-E(Vcon);
%     E(E>20)=20; E(E<-20)=-20;
             
    op=sum(E(Vcon),2)+V;
    op(op>20)=20; op(op<-20)=-20;

% dec=sign(sum(E(Vcon),2)+V);
end