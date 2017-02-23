% clc
clear all
global dv dc w L M
dv=3;
dc=10;
w=3; %Works only for w=3 and dv=3
M=50;
L=16;
U=2; % Number of users.

llr_max=20; 
maxIters=100;
maxSims=1000;
chan='noiseless';
sigVec=0;

%% For two different graphs
[Vcon1,Ccon1,eMax]=ECon_SC();
[Vcon2,Ccon2,eMax]=ECon_SC();
N=length(Vcon1(:,1));
m=length(Ccon1(:,1)); %Number of check nodes

%% 
fprintf('Sys. Params are N=%d,m=%d,dc=%d,dv=%d, rate=%4.3f via M=%d\n',N,m,dc,dv,1-m/N,M);

for sigIter=1:length(sigVec)
    sig=sigVec(sigIter);
    fprintf('Channel is %s and sig=%f\n',chan,sig');

    num_err=0;  sim_cnt=0; blockErr=0;
    tic

 while blockErr<10 && sim_cnt<maxSims
    
    c=rand(2,N)<0.5;
    E=zeros(1,eMax+1);
    E(Vcon1)=repmat(c(1,:)',1,dv); P(1,:)=mod(sum(E(Ccon1),2),2);
    E=zeros(1,eMax+1);
    E(Vcon2)=repmat(c(2,:)',1,dv); P(2,:)=mod(sum(E(Ccon2),2),2);
    
    x=1-2*c;
    y=x(1,:)+x(2,:)+sig*randn(1,N);

    xhat=BP_two_user_MAC_2(y,chan,[Vcon1;Vcon2],[Ccon1;Ccon2],eMax,maxIters,P',sig);

    num_err=numel(find(xhat-x))+num_err;
    blockErr=~isempty(find(xhat-x,1))+blockErr;
    sim_cnt=sim_cnt+1;

    if numel(find(xhat-x))==1
       fprintf('For debugging') 
    end
    
    if mod(sim_cnt,50)==0
      fprintf('#err bits= %d/%d Sims. %3.2f sec\n',num_err,sim_cnt,toc);
    end
  end
    
fprintf('# bits in error= %d after %d Simulations\n',num_err,sim_cnt);
errMat(sigIter,1:3)=[num_err blockErr sim_cnt];
end