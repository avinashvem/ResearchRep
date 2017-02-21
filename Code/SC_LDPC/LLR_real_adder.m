function [llr1]=LLR_real_adder(y,llr_max,channel,U,sig)
N=numel(y);
if U==2
   if strcmp(channel,'noiseless')
       llr1=zeros(numel(y),1);        
       pos= y==2; neg= y==-2;
       llr1(pos)=llr_max; llr1(neg)=-llr_max; 
    elseif strcmp(channel,'Gaussian')
       llr1=4*y/sig^2;
   end
end




end