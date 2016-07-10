function [Vout]=GenProd(Vin,beta)

n=size(Vin,2); m=size(Vin,1);
Vin=tanh(beta*Vin);

Vf=cumprod(Vin(:,1:(n-1)),2); Vf=[ones(m,1) Vf];
Vb=cumprod(Vin(:,n:-1:2),2); Vb=Vb(:,n-1:-1:1); Vb=[Vb ones(m,1)];
Vout=Vf.*Vb;

assert(all(size(Vin)==size(Vout)));