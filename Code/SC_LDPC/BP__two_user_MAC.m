function op=BP__two_user_MAC(Vchan,Vcon,Ccon,eMax,Parity,shiftPattern)
% dv and dc stand for Max Variable and Check node degrees respectively .This is intended for a Irregular degree profile too
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

if numel(Vchan(1,:))~=2
   Vchan=Vchan'; 
end

E1=zeros(e,1);
E1(e+1)=Inf;
E2=E1;

Ccon(Ccon>e)=e+1;

E_prev1=zeros(length(E1),1);
E_prev2=zeros(length(E2),1);

E1(Vcon)=repmat(Vchan(:,1),1,dv);
E2(Vcon)=repmat(Vchan(:,2),1,dv);
%% Message Passing
% for i=1:iter
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
% end
% dec=sign(sum(E(Vcon),2)+V);
end