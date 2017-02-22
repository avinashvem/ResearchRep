function [llr_out]=LLR_real_adder(y,llr_in,channel,U,sig)
global llr_max
l1=llr_in(:,1);
l2=llr_in(:,2);
if U==2
   if strcmp(channel,'noiseless')
%        llr_out=zeros(numel(y),1);        
%        pos= y==2; neg= y==-2;
%        llr_out(pos)=llr_max; llr_out(neg)=-llr_max; 
       llr_out(:,1)=-l2.*sign(2-abs(y))+sign(y).*llr_max;
       llr_out(:,2)=-l1.*sign(2-abs(y))+sign(y).*llr_max;

    elseif strcmp(channel,'Gaussian')
       llr_out(:,1)=log( (1+exp(l2+2*(y-1)/sig^2))./(exp(l2)+exp(-2*(y+1)/sig^2)) );
       llr_out(:,2)=log( (1+exp(l1+2*(y-1)/sig^2))./(exp(l1)+exp(-2*(y+1)/sig^2)) );
       llr_out(llr_out>20)=20;  llr_out(llr_out<-20)=-20;
   end
end
end