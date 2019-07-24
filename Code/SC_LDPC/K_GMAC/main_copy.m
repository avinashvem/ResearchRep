% clear all
addpath('../PEG/')
global M repeat_flag repeat interl deinterl

Ka=50;
N=30000;
M=400;

L=[0.4 0.2 0.4;2 3 4];R=[1;4];
base_rate=1/4;

EsN0_vec=-3.75;
maxIters=500;

lavg=sum(L(1,:).*L(2,:));
ravg=sum(R(1,:).*R(2,:));
assert(ravg*(1-base_rate)==lavg);
repeat_flag=0;
repeat=1;

%--Michael's PEG construction------------------
% [Vcon,Ccon,eMax]=LDPC_PEG_cnstr(M,round(M*(1-base_rate)),L,R,repeat_flag);

fprintf('L(x)=');
for i=1:numel(L(1,:))
   fprintf('%3.2f x^%d+',L(1,i),L(2,i));
end
fprintf('\n R(x)=');
for i=1:numel(R(1,:))
   fprintf('%3.2f x^%d+',R(1,i),R(2,i));
end
fprintf('\n')

n=length(Vcon(:,1)); dv=size(Vcon,2);
m=length(Ccon(:,1)); %Number of check nodes

interl=zeros(N,Ka);

EbN0_vec=EsN0_vec+10*log10(n/(n-m));
sigVec=sqrt(0.5./10.^(EsN0_vec./10));
llr_max=20; 
maxSims=20;

%-------------Simulations-----------------------------------------
fprintf('Sys.params: N=%d, (k,n)=(%d,%d), Ka=%d \n',N,n-m,n,Ka);
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
    while blockErr<10 && sim_cnt<maxSims
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
    
    fprintf('#err bits/blocks= %d/%d, Berr=%4e in %d Simulations\n',num_err,blockErr,blockErr/sim_cnt,sim_cnt);
%   fprintf('#err bits= %d,in %d Simulations\n',num_err,sim_cnt);
%   errMat(sigIter,1:3)=[num_err bitErr sim_cnt];
end

%-----Simulation Results--------------------