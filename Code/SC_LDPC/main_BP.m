% clc
clear all
global dv dc w L M
dv=3;
dc=11;
w=3; %Works only for w=3 and dv=3
M=50;
L=16;
U=2; % Number of users.
shift=round(4.5*M);

llr_max=20; 
maxIters=100;
maxSims=1000;

[Vcon,Ccon,eMax]=ECon_SC();
N=length(Vcon(:,1));
m=length(Ccon(:,1)); %Number of check nodes

%% 
fprintf('Sys. Params are N=%d,m=%d,dc=%d,dv=%d, rate=%4.3f via M=%d\n',N,m,dc,dv,1-m/N,M);

E=zeros(1,eMax+1);
shiftPattern=[shift+1:N 1:shift];
inv_shiftPattern=[N-shift+1:N 1:N-shift];

chan='Gaussian';
sigVec=[0.34];
for sigIter=1:length(sigVec)
sig=sigVec(sigIter);
fprintf('Channel is %s and sig=%f\n',chan,sig');

num_err=0;  sim_cnt=0; blockErr=0;
tic
while blockErr<2 && sim_cnt<maxSims
for u=1:U
    c(u,:)=rand(1,N)<0.5;
    E=zeros(1,eMax+1);
    E(Vcon)=repmat(c(u,:)',1,dv);
    P(u,:)=mod(sum(E(Ccon),2),2);
end
x=1-2*c;
y=x(1,:)+x(2,shiftPattern)+sig*randn(1,N);

xhat=BP_two_user_MAC(y,chan,Vcon,Ccon,eMax,maxIters,P',shiftPattern',sig);

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

%% ----Simulation Results-----
% dv=3, dc=9,L=16,w=3
% M    (N,m)      (bErr,Berr) #Sims  sig
% 50 (2400,900)     (0,0)      500   0

% dv=3, dc=10,L=16,w=3, rate=0.6625
% M    (N,m)      (bErr,Berr) #Sims  sig
% 50 (8000,2700)    (0,0)      500    0
% 50 (8000,2700)   (6.1e4,20)   20   0.50
% 50 (8000,2700)   (1.7e4,20)   52   0.40
% 50 (8000,2700)   (1.5e4,20)   98   0.38
% 50 (8000,2700)   (1.2e4,20)  235   0.36
% 50 (8000,2700)   (1.3e4,16)  500   0.34
% 50 (8000,2700)   (2.5e3,4)   500   0.32
% 50 (8000,2700)   (1.1e4,8)  1500   0.30
% 50 (8000,2700)   (1.4e3,01)  1e3   0.29
%500 (8e4, 27e3)    (0,0)       250   0.34

% dv=3, dc=11,L=16,w=3, rate=0.6932
% M    (N,m)      (bErr,Berr) #Sims  sig
% 50 (8800,2700)    (7e4,20)    20    0

% dv=4, dc=14,L=16,w=3, rate=0.679
% M    (N,m)      (bErr,Berr) #Sims  sig
%90 (1e4,3240)    (5.5e4,20)   20    0