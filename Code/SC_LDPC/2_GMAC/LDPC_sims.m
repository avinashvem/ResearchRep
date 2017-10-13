% clc
% clear all
global dv dc M repeat_flag repeat
dv=3;
dc=6;
M=2000;
L=[1;3];R=[1;6];
% L=[0.5 0.5;2 3];  R=[1;5];%0.5 0.5;5 6];
base_rate=1/2;

lavg=sum(L(1,:).*L(2,:));
ravg=sum(R(1,:).*R(2,:));
assert(ravg*base_rate==lavg);
repeat_flag=1;
repeat=3;


%--Michael's PEG construction------------------
if repeat_flag
    [Vcon,Ccon,eMax]=LDPC_PEG_cnstr(M,round(M*(1-base_rate)),L,R,repeat_flag);
else
    [Vcon,Ccon,eMax]=LDPC_PEG_cnstr(M,round(M*(1-base_rate)),[1;dv],repeat_flag);
end
% [Vcon,Ccon,eMax]=ECon_LDPC(M,round(M*(1-base_rate)),L,R,repeat_flag);
% [Vcon,Ccon,eMax]=IRA_rate_half(M,repeat_flag);

N=length(Vcon(:,1)); dv=size(Vcon,2);
m=length(Ccon(:,1)); %Number of check nodes

U=2; % Number of users.
EsN0_vec=[-6.25 -5.85 -5.7];%6.20 -6.0 -5.80 -5.5 -5.25 -5.0 -4.9];

sigVec=sqrt(0.5./10.^(EsN0_vec./10));
llr_max=20; 
maxIters=200;
maxSims=10000;

%-------------Simulations-----------------------------------------
fprintf('Sys. Params are N=%d,m=%d,M=%d rate=%4.3f\n',N,m,M,1-m/N);
fprintf('L(x)=');
for i=1:length(L(1,:))
    fprintf('%3.2f x^%d +',L(1,i),L(2,i));
end
fprintf('... R(x)=');
for i=1:length(R(1,:))
    fprintf('%3.2f x^%d +',R(1,i),R(2,i));
end
fprintf('\n');

shiftPattern=randperm(N);
chan='Gaussian';

for sigIter=1:length(sigVec)
    sig=sigVec(sigIter);
    fprintf('Channel is %s and sig=%f or EsNo=%3.2f dB\n',chan,sig',10*log10(1/2/sig^2));

    num_err=0;  sim_cnt=0; blockErr=0;
    tic
    while blockErr<10 && sim_cnt<maxSims
        for u=1:U
            if repeat_flag
                c(u,1:M)=rand(1,M)<0.5;
                for rpt=1:repeat-1
                    c(u,rpt*M+1:(rpt+1)*M)=c(u,1:M);
                end
            else
                c(u,:)=rand(1,N)<0.5;
            end
            E=zeros(1,eMax+1);
            E(Vcon)=repmat(c(u,:)',1,max(dv));
            E(eMax+1)=0;
            P(u,:)=mod(sum(E(Ccon),2),2);
        end
        
        x=1-2*c;
        y=x(1,:)+x(2,shiftPattern)+sig*randn(1,N);
        
        xhat=BP_2GMAC(y,chan,Vcon,Ccon,eMax,maxIters,P',shiftPattern',sig);
        
%         if numel(find(xhat-x))>10
           num_err=num_err+numel(find(xhat-x));
           blockErr=~isempty(find(xhat-x,1))+blockErr;
%         end
        
        
        sim_cnt=sim_cnt+1;

        if numel(find(xhat-x))==1
           fprintf('For debugging') 
        end
    
        if mod(sim_cnt,500)==0
          fprintf('#err bits/blocks= %d/%d in %d Sims. %3.2f sec\n',num_err,blockErr,sim_cnt,toc);
        end
    end
    
    fprintf('#err bits/blocks= %d/%d, Berr=%4e inr %d Simulations\n',num_err,blockErr,blockErr/sim_cnt,sim_cnt);
    errMat(sigIter,1:3)=[num_err blockErr sim_cnt];
end

%-----Simulation Results--------------------
%{
% (3,6) repeated 3 times. TotalRate=1/6
% Es_N0 Berrate  Berr Sims   N   berr
% -4.0  4.90e-2   10   204  600  106  Berrate=1.44e-3 with 2-cycles fixed
% -3.5  1.20e-3    6  5000  600  130  Avg 22bits in Error/Berr
% -3.0  1.36e-3   10  7350  600  105

% -5.0  1.56e-2    10  640 1200  868
% -4.7  3.30e-3    10 3026 1200  292
% -4.4  3.155e-3   10 3169 1200  138
% -4.2  3.24e-3    10 3084 1200  113
% -4.0  3.00e-3    10 3287 1200  131
%-3.85  8.33e-4     5 6000 1200   69 
%-3.7   0           0 2500 1200   0

%-5.9   1.15e-1    10   87 3000  2611
%-5.7   5.02e-2    10  199 3000  2651
%-5.5   2.35e-2    10  424 3000  2651
%-5.3   6.66e-4     2 3000 3000    45
%-5.0   0          0  2500 3000   0

% (3,6) repeated 2 times. TotalRate=1/4
% Es_N0 Berrate  Berr Sims   N   berr
% -3.0  2.79e-2   10   358  400  710 
% -2.5  7.13e-3   10  1402  400  473
% -2.3  2.30e-3   10  4348  400  404
% -2.0  1.44e-3   10  6925  400  165
% -1.7  1.07e-3   10  9317  400  138

%After Fixing 2,4-cycles
% (3,6) repeated 2 times. TotalRate=1/4
% Es_N0 Berrate  Berr Sims   N   berr
% -1.7  7.00e-4   7  10000  400   79
% -1.5  2.00e-4   2  10000  400   47

% cap_Sh(1/6,T=2)=-8.32dB
% inv_random_codinng_GMAC results for P1=(0.11:0.04:0.99)*P
%   N   K Berrate  Es_N0 
%  600 100 1e-4    -6.48
%  600 100 1e-3    -6.9922
% 1200 200 1e-4    -7.6562
% 1200 200 1e-3    -8.0477
% 3000 500 1e-3    -8.8281

% inv_random_codinng_GMAC results for P1=(0.11:0.04:0.99)*P
%  N   K Berrate   Es_N0 
% 400 100 1e-4    -4.4531
% 400 100 1e-3    -4.9250
%}