function [Vcon,Ccon,e]=ECon_SC()
global dv dc w L M
% dv=3;
% dc=9;
% w=3; %Works only for w=3 and dv=3% 
% M=3000;
% L=16;

 m=3; % m = no of check nodes in base protograph at each position
 n=10; %It should be m*dc=n*dv


Vposn=n*M;
Cposn=m*M; %No. of check nodes at each position
ePosn=n*M*dv;

Vtot=Vposn*L;   
Ctot=Cposn*(L+w-1);
e=Vtot*dv;%alternately e1*L
e_dummy=Ctot*dc;

Vcon=reshape(1:e,dv,Vtot)';

V=(e+1)*ones(w,ePosn);
V(1,:)=randperm(ePosn);
e2=ePosn/w;
if round(e2)~=e2
   fprintf('Error.Number of edges not divisible by w.')
end

for i=1:L+w-1
    temp_e=[V(1,1:e2) V(2,e2+1:2*e2) V(3,2*e2+1:3*e2)];
    Ccon((i-1)*Cposn+1:i*Cposn,:)=reshape(temp_e,[],dc);
    V(2:3,:)=V(1:2,:);
    V(1,:)=i*ePosn+randperm(ePosn);
end
    
Ccon(Ccon>(e+1))=e+1;