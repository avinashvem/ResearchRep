% clc
clear all
global dv dc w L M
dv=3;
dc=10;
w=3; %Works only for w=3 and dv=3
M=50;
L=16;
U=2; % Number of users.
shift=4.5*M;

llr_max=20; 
maxIter=500;
maxSims=1000;
[Vcon,Ccon,eMax]=ECon_SC();
N=length(Vcon(:,1));
m=length(Ccon(:,1)); %Number of check nodes

chan='noiseless';
sig=0;
fprintf('Sys. Params are N=%d,m=%d,dc=%d,dv=%d\n',N,m,dc,dv);

E=zeros(1,eMax+1);
shiftPattern=[shift+1:N 1:shift];
inv_shiftPattern=[N-shift+1:N 1:N-shift];

num_err=0;  sim_cnt=0; blockErr=0;
tic
while num_err<500 && sim_cnt<maxSims
for u=1:U
    c(u,:)=rand(1,N)<0.5;
    E=zeros(1,eMax+1);
    E(Vcon)=repmat(c(u,:)',1,dv);
    P(u,:)=mod(sum(E(Ccon),2),2);
end
x=1-2*c;
y=x(1,:)+x(2,shiftPattern)+sig*randn(1,N);

V1Chan=LLR_real_adder(y,llr_max,chan,U,sig);
V2Chan=V1Chan(inv_shiftPattern);


V1_prev=V1Chan;  V2_prev=V2Chan;
zeroIdx= find(V1Chan==0);

% BP__two_user_MAC([V1Chan V2Chan],Vcon,Ccon,eMax,Parity,shiftPattern)

for i=1:maxIter    
    V1=BP_Iregular_BMSC(V1_prev,Vcon,Ccon,eMax,P(1,:)');
    V2=BP_Iregular_BMSC(V2_prev,Vcon,Ccon,eMax,P(2,:)');

    if strcmp(chan,'noiseless')
        V1_prev=V1Chan;   V2_prev=V2Chan;
        V2_shifted=V2(shiftPattern);
    
        V1_prev(zeroIdx)=-V2(shiftPattern(zeroIdx));   
        V2_prev(shiftPattern(zeroIdx))=-V1(zeroIdx);
    end
   
    if max(V1_prev-V1)<0.001 && max(V2_prev-V2)<0.001
        xhat=sign([V1+V1_prev V2+V2_prev])';
        break 
    end   
    %{
    if mod(i,100)==0
        xhat=sign([V1+V1_prev V2+V2_prev])';
        fprintf('Number of errors is %d\n', numel(find(xhat-x)));
    end
    %}
end

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


%% ----Simulation Results-----
% dv=3, dc=9,L=16,w=3
% M    (N,m)      (bErr,Berr)   #Sims
% 50 (2400,900)     (0,0)      1000
%100 (4800,1800)    (6,5)        1000