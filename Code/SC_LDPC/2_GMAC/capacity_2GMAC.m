function cap=capacity_2GMAC(ip,rateFlag)
%ip= sig or the rate. If rateFlag==1 then ip is the single user rate.
%if rateFlag==0 then the ip is the sigma(std_dev of Gaussian noise).
%For a given rate/sig, compute the capacity of 2 user real adder Gaussian MAC channel. 
%y=x1+x2+n where x1,x2 in {-1,+1} and n is AWGN(0,sig^2).


if rateFlag==0
    sig=ip;
    f1 = @(x) exp(-((x+2).^2/2/sig^2));
    f0 = @(x) exp(-((x.^2)/2/sig^2));
    f3 = @(x) exp(-((x-2).^2/2/sig^2));
    f=  @(x) (0.25/sqrt(2*pi*sig^2))*(f1(x)+2*f0(x)+f3(x));
    g=@(x) -f(x).*log2( 0.25*(f1(x)+2*f0(x)+f3(x)) );

    cap=integral(g,-25*sig,25*sig)-0.5*log2(exp(1));

else
    rate=ip;
    % function EsN0=generate_capacity(rate)
    % For a given single user rate, compute the capacity int terms of min EsN0(dB) required for a 
    % 2 user real adder Gaussian MAC channel. 
    %y=x1+x2+n where x1,x2 in {-1,+1} and n is AWGN(0,sig^2).

    sig_high=2;
    sig_low=0;

    while sig_high-sig_low>0.001
        sig=(sig_high+sig_low)/2;
        if capacity_2user_GMAC(sig,0)>2*rate
            sig_low=sig;
        else
            sig_high=sig;
        end
    end
    EsN0=10*log10(1/2/sig^2);
    fprintf('For a desired sum-rate of %3.2f, min Es/N0 is %3.2f dB\n',2*rate,EsN0);
end
end

%{
sigVec=0.1:0.01:0.9;
SNR=10*log10(1./(2*(sigVec.^2)));
cap=zeros(size(sigVec));

for sigIter=1:length(sigVec)
    cap(sigIter)=capacity_2user_GMAC(sigVec(sigIter));
end
plot(SNR,cap,'b')
%}