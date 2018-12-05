% clear all
addpath('../PEG/')
global M repeat_flag repeat interl deinterl

Ka=25;
N=30000;
M=200;
L=[0.5 0.5;2 3];R=[1;5];
base_rate=1/2;

EsN0_vec=-1:0.2:1;
maxIters=500;
maxSims=20;

lavg=sum(L(1,:).*L(2,:));
ravg=sum(R(1,:).*R(2,:));
assert(ravg*(1-base_rate)==lavg);
repeat_flag=0;
repeat=1;

%--Michael's PEG construction------------------
[Vcon,Ccon,eMax]=LDPC_PEG_cnstr(M,round(M*(1-base_rate)),L,R,repeat_flag);

fprintf('L(x)=');
for i=1:numel(L(1,:))
   fprintf('%3.2f x^%d+',L(1,i),L(2,i));
end
fprintf('\n R(x)=');
for i=1:numel(R(1,:))
   fprintf('%3.2f x^%d+',R(1,i),R(2,i));
end
fprintf('\n');

n=length(Vcon(:,1)); dv=size(Vcon,2);
m=length(Ccon(:,1)); %Number of check nodes
interl=zeros(N,Ka);

EbN0_vec=EsN0_vec+10*log10(n/(n-m));
sigVec=sqrt(0.5./10.^(EsN0_vec./10));

llr_max=20;  
%-------------Simulations-----------------------------------------
fprintf('Sys.params: N=%d,Ka=%d, (k,n)=(%d,%d) \n',N,Ka,n-m,n);

for i=1:Ka
    interl(:,i)=reshape(randperm(N),[],1);
    deinterl(interl(:,i),i)=1:N;
end
chan='Gaussian';

for sigIter=1:length(sigVec)
    sig=sigVec(sigIter);
    fprintf('%s ch. and sig=%f, Es/N0=%3.2f dB,Eb/N0=%3.2f\n',chan,sig',EsN0_vec(sigIter),EbN0_vec(sigIter));

    num_err=0;  sim_cnt=0; blockErr=0;
    tic
    while num_err<20 && sim_cnt<maxSims
        c=zeros(N,Ka);
        x=c;
        for u=1:Ka
            cw(1:n,u)=(rand(n,1)<0.5);
            c(1:n,u)=1-2*cw(1:n,u);
            x(:,u)=c(interl(:,u),u);
            
            E=zeros(1,eMax+1);
            E(Vcon)=repmat(cw(1:n,u),1,max(dv));
            E(eMax+1)=0;
            P(:,u)=mod(sum(E(Ccon),2),2);
        end
        
    y=sum(x,2)+sig*randn(N,1);
    xhat=BP_KGMAC(y,Vcon,Ccon,eMax,maxIters,P,sig);
        
%------Book-keeping for Errors--------
       for u=1:Ka
        num_err=num_err+(~isempty(find(xhat(:,u)-x(:,u),1)));
       end
       
       sim_cnt=sim_cnt+1;
    
       if mod(sim_cnt,1)==0
           fprintf('#err users/total= %d/%d---%3.2f sec\n',num_err,Ka*sim_cnt,toc);
       end
    end
    
    fprintf('#err bits= %d, berr=%4e in %d Simulations\n',num_err,num_err/Ka/sim_cnt,sim_cnt);
%     fprintf('#err bits= %d,in %d Simulations\n',num_err,sim_cnt);
%     errMat(sigIter,1:3)=[num_err bitErr sim_cnt];
end

%-----Simulation Results--------------------

%-----------N=10,000------------------------------
% Ka=50;L=x^3;R=x^6;k=100;n=200;
% Eb/N0 Errs Sims  Es/N0  eps
% 3.51   1    15    0.50
% 2.86   0     2   -0.15  0.00
% 2.81   92    4   -0.20  0.4600
% 2.76   848  20   -0.25  0.8480  
% 2.51   489  10   -0.50  0.9780
% 2.01   999  20   -1.00  0.9999

% Ka=50; L=0.5 x^2+0.5x^3;R=x^5;k=100;n=200;
% Eb/N0 Errs Sims Es/N0  eps
% 3.01   9    10   0.00 1.8e-02
% 2.56  18     6  -0.45 6.0e-02

%------Ka=50; N=30,000-----------

% L=x^3;R=x^6;k=100;n=200;
% Eb/N0 Errs Sims Es/N0  eps
% 3.01   0    2   0.00   0.00
% 2.26   16   5  -0.75   6.40e-2

% L=0.5x^2+0.5x^3;R=x^5;k=100;n=200;
% Eb/N0 Errs Sims Es/N0  eps
% 2.06   49    8  -0.95  1.23e-01
% 2.11   30    8  -0.90  7.50e-02

% Rate-1/4 -> k=100;n=400;
% L=4x^2+2x^3+4x^4; R=x^4;
% Eb/N0 Errs Sims Es/N0  eps
% 4.02   0    3   0.00  0.00
% 2.92   0    4  -3.10  0.00
% 2.27   3    3  -3.75  2.00e-02
% 2.02   9    4  -4.00  4.50e-02


%------Ka=25; N=30,000-----------

% L=x^3;R=x^6;k=100;n=200;
% Eb/N0 Errs Sims Es/N0  eps
% 

% L=0.5x^2+0.5x^3;R=x^5;k=100;n=200;
% Eb/N0 Errs Sims Es/N0   eps
% 1.81   49    10  -1.20  1.96e-01 
% 1.91   37    10  -1.10  1.48e-01 
% 2.01   21    5   -1.00  1.68e-01 
% 2.21   20    8   -0.80  1.00e-01 
% 2.41   20    18  -0.60  4.44e-02
% 2.61   14    20  -0.40  2.80e-02
% 2.81    8    20  -0.20  1.60e-02
% 3.01    4    20   0.0   8.00e-03



% Rate-1/4 -> k=100;n=400;
% L=4x^2+2x^3+4x^4; R=x^4;
% Eb/N0 Errs Sims Es/N0  eps
%