function cap=capacity_Tuser_GMAC(ip,T,rateFlag)
%ip= sig or the rate. If rateFlag==1 then ip is the single user rate.
%if rateFlag==0 then the ip is the sigma(std_dev of Gaussian noise).
%For a given rate/sig, compute the capacity of 2 user real adder Gaussian MAC channel. 
%y=x1+x2+n where x1,x2 in {-1,+1} and n is AWGN(0,sig^2).

if mod(T,2)~=0
   error('Input even values of T')
end
   
if rateFlag==0
    sig=ip;
    f=@(x) 0;
    for i=0:T/2
       f = @(x) f(x)+nchoosek(T,T/2-i)*exp(-((x-2*i).^2/2/sig^2));
    end
    for i=-T/2:-1
       f = @(x) f(x)+nchoosek(T,T/2+i)*exp(-((x-2*i).^2/2/sig^2));
    end
    
    g=@(x) -2^(-T)/sqrt(2*pi*sig^2)*f(x).*log2( 2^(-T)*f(x));
    sum_rate=integral(g,-25*sig,25*sig)-0.5*log2(exp(1));
    cap=sum_rate/T;
else
    rate=ip;
    % function EsN0=generate_capacity(rate)
    % For a given single user rate, compute the capacity int terms of min EsN0(dB) required for a 
    % 2 user real adder Gaussian MAC channel. 
    %y=x1+x2+n where x1,x2 in {-1,+1} and n is AWGN(0,sig^2).

    sig_high=3;
    sig_low=0;

    while sig_high-sig_low>0.001
        sig=(sig_high+sig_low)/2;
        if capacity_Tuser_GMAC(sig,0)>rate
            sig_low=sig;
        else
            sig_high=sig;
        end
    end
    EsN0=10*log10(1/2/sig^2);
    fprintf('For a desired sum-rate of %3.2f, min Es/N0 is %3.2f dB\n',2*rate,EsN0);
end
end