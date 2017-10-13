% clc
addpath('../')
clear all
global dv dc w L M
dv=4;
dc=6;
w=3; %Works only for w=3 and dv=3
M=100;
L=16;
U=2; % Number of users.
shift=round(4.5*M);
EsN0_vec=[-4.500 -4.300 -3.90 -3.80 -3.70 -3.60];
sigVec=sqrt(0.5*(10.^(-EsN0_vec./10)));

llr_max=20; 
maxIters=500;
maxSims=1000;

[Vcon,Ccon,eMax]=ECon_SC();
N=length(Vcon(:,1));
m=length(Ccon(:,1)); %Number of check nodes

%----
fprintf('Sys. Params are N=%d,m=%d,dc=%d,dv=%d, rate=%4.3f via M=%d\n',N,m,dc,dv,1-m/N,M);

E=zeros(1,eMax+1);
% shiftPattern=[shift+1:N 1:shift];
% inv_shiftPattern=[N-shift+1:N 1:N-shift];
shiftPattern=randperm(N);
inv_shiftPattern(shiftPattern)=1:N;


chan='Gaussian';

for sigIter=1:length(sigVec)
sig=sigVec(sigIter);
fprintf('Channel is %s and sig=%f or EsNo=%3.2f dB\n',chan,sig',10*log10(1/2/sig^2));

num_err=0;  sim_cnt=0; blockErr=0;
tic
while blockErr<10 && sim_cnt<maxSims
for u=1:U
    c(u,:)=rand(1,N)<0.5;
    E=zeros(1,eMax+1);
    E(Vcon)=repmat(c(u,:)',1,dv);
    P(u,:)=mod(sum(E(Ccon),2),2);
end
x=1-2*c;
y=x(1,:)+x(2,shiftPattern)+sig*randn(1,N);

xhat=BP_2GMAC(y,chan,Vcon,Ccon,eMax,maxIters,P',shiftPattern',sig);

num_err=numel(find(xhat-x))+num_err;
blockErr=~isempty(find(xhat-x,1))+blockErr;
sim_cnt=sim_cnt+1;

if numel(find(xhat-x))==1
   fprintf('For debugging') 
end
    
if mod(sim_cnt,50)==0
  fprintf('#err bits/blocks= %d/%d in %d Sims. %3.2f sec\n',num_err,blockErr,sim_cnt,toc);
end
end
    
fprintf('#err bits/blocks= %d/%d after %d Simulations\n',num_err,blockErr,sim_cnt);
errMat(sigIter,1:3)=[num_err blockErr sim_cnt];
end
%------Simulation Results------------------------------------------------------------------
%{
% (6,8) dv=6, dc=8,L=16,w=3, rate=0.156
%  M     (N,m)    (bErr,Berr)  #Sims Es_N0  sig
% 100 (6400,5400)    (0,0)      1000  -3.0  0.9988
% 10   (640,540)    (774,10)     506  -0.0  0.7071
% 10   (640,540)    (376,5)     1000  1.00  0.6302
% 10   (640,540)    (75,1)      1000  1.75  0.6302


% (3,4) dv=3, dc=4,L=16,w=3, rate=0.156
% M=100; (N,m)=(64,54)*1e2
% (bErr,Berr) #Sims Es_N0  sig
%   (18,1)     1000  --    1.20
%   (0,0)      1000  -3.0  0.9988

% (3,6); (N,m)=(8000,4500),M=250 L=16,w=3, rate=0.4375
%     (bErr,Berr)  #Sims    sig
%     (4.3e4,50)    196    0.750
%     (4.6e4,50)    239    0.747
%     (4.3e4,50)    404    0.744
%     (3.9e4,50)    481    0.741
%     (3.0e4,50)    741    0.738
%     (2.5e4,39)   1000    0.735
%     (1.5e4,24)   1000    0.732
%      (5e3,10)     782    0.729
%      (5e3,6)     1000    0.726
%      (2e3,3)     1000    0.723
%      (294,2)      38     0.73
%       (0,0)      1000    0.72
%       (0,0)      1000    0.71

% (3,6); (N,m)=(32000,1800),M=1000 L=16,w=3, rate=0.4375
% (bErr,Berr) #Sims  sig
% (3.2e5,50)     55   0.790
% (2.8e5,50)    171   0.780
% (2.3e5,45)   1000   0.770
% (1.0e5,18)   1000   0.767
% (4.0e4,8)    1000   0.764
% (2.0e4,4)    1000   0.762
% (7006,1)     1000   0.760
%   (0,0)      1000   0.756

% (3,6); (N,m)=(96000,54000),M=3000 L=16,w=3, rate=0.4375
%   (bErr,Berr) #Sims  sig
%     (N,2)        2   0.80
%   (1.8e5,10)    12   0.79
%   (1.5e5,10)    88   0.78
%   (1.8e4,1)    650   0.77
%    (1178,2)    411   0.76
%      (0,0)    1050   0.75
%      (0,0)     600   0.70

% (3,9); (N,m)=(7200,2700),M=150 L=16,w=3, rate=0.625
% (bErr,Berr) #Sims  sig
% (4.7e4,40)    107   0.50
% (4.7e4,40)    199   0.49
% (4.7e4,40)    674   0.48
% (4.7e4,50)    788   0.48
% (1.0e4,22)   1000   0.47
% (6.2e3,12)   1000   0.46
% (4.7e3,8)    1000    0.45
% (2.3e3,4)    1000    0.44
% (1080,2)     1000    0.43
% (1614,2)     1000    0.42

% (3,9); (N,m)=(30240,11340),M=630 L=16,w=3, rate=0.625
% (bErr,Berr)  #Sims    sig
%  (7.0e4,10)     10   0.530
%  (6.2e4,10)     22   0.520
%  (1.8e5,40)    790   0.510
%  (9.6e4,20)    596   0.508
%  (6.6e4,13)    750   0.505
%  (2.2e4,5)     750   0.503
%  (8700,2)     1000   0.500
%    (0,0)      1000   0.49
 
% (3,9); (N,m)=(33600,12600),M=700 L=16,w=3, rate=0.625 
% (bErr,Berr) #Sims  sig
%  (265,2)      2   0.500
% (3.7e2,2)     2   0.490
% (1400,2)    473   0.480
% (1241,1)    1e3   0.477
%  (0,0)      1e3   0.470

% (3,9); (N,m)=(96000,36000),M=2000 L=16,w=3, rate=0.625 
%(bErr,Berr) #Sims  sig
%  (8e3,2)      2   0.500
% (1249,2)     92   0.490
%  (0,0)      1e3   0.480
% (199,2)     765   0.486
%  (0,0)      1e3   0.483

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
%}

