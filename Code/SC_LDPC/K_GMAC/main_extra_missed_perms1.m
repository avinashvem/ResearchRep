%{
The extension to main_extra_missed_perms is that one simulation is
performed for all input Es_N01 values, and input list of permutations and
decoded list of permutations.
%}
function errMat=main_extra_missed_perms1(chcode,params,Xcs,Xcshat)
global M repeat_flag repeat

N=chcode.N;
M=chcode.M;
L=chcode.L; R=chcode.R;
base_rate=chcode.base_rate;
repeat_flag=chcode.repeat_flag;
repeat=chcode.repeat;

Vcon=chcode.Vcon;
Ccon=chcode.Ccon;
eMax=chcode.eMax;

Ka=params.Ka;
EsN0_vec=params.EsN0_vec;


maxIters=500;

lavg=sum(L(1,:).*L(2,:));
ravg=sum(R(1,:).*R(2,:));
assert(ravg*(1-base_rate)==lavg);

n=length(Vcon(:,1)); dv=size(Vcon,2);
m=length(Ccon(:,1)); %Number of check nodes


total_indices=union(Xcs,Xcshat);
Ktotal=numel(total_indices);
interl_map=sort(total_indices);

Xcs=reshape(Xcs,[],1);
tmp_Xcshat=Xcshat;
unfound_idx=[];
for i=1:Ka
    if ismember(Xcs(i),tmp_Xcshat)
        idx=find(tmp_Xcshat==Xcs(i),1);
        tmp_Xcshat(idx)=0;
        Xcshat(i)=Xcs(i);    
    else
        unfound_idx(end+1)=i;
    end
end
unfound_idx=[unfound_idx Ka+1:length(Xcshat)];
iter=1;
tmp_Xcshat=tmp_Xcshat(find(tmp_Xcshat));
for i=1:length(unfound_idx)
   Xcshat(unfound_idx(i))=tmp_Xcshat(iter);
   iter=iter+1;
end
    
    
Xcs_int=zeros(length(Xcs),1);
Xcshat_int=zeros(length(Xcshat),1);
for i=1:numel(Xcs)
    Xcs_int(i)=find(interl_map==Xcs(i));
end

for i=1:numel(Xcshat)
    Xcshat_int(i)=find(interl_map==Xcshat(i));
end

interl=zeros(N,Ktotal);
deinterl=interl;

EbN0_vec=EsN0_vec+10*log10(n/(n-m));
sigVec=sqrt(0.5./10.^(EsN0_vec./10));
errMat=-1*ones(length(EbN0_vec),4);
%-------------Simulations-----------------------------------------
for i=1:Ktotal
    interl(:,i)=reshape(randperm(N),[],1);
    deinterl(interl(:,i),i)=1:N;
end

for sigIter=1:length(sigVec)
    sig=sigVec(sigIter);
    

    num_err=0;  sim_cnt=0;
    
    c=zeros(N,Ka);
    x=c;
    for u=1:Ka
        cw(1:n,u)=(rand(n,1)<0.5);
        c(1:n,u)=1-2*cw(1:n,u);
        x(:,u)=c(interl(:,Xcs_int(u)),u);

        E=zeros(1,eMax+1);
        E(Vcon)=repmat(cw(1:n,u),1,max(dv));
        E(eMax+1)=0;
        P(:,u)=mod(sum(E(Ccon),2),2);
    end
    
     
    y=sum(x,2)+sig*randn(N,1);
    xhat=BP_KGMAC_extraperms(y,Vcon,Ccon,eMax,maxIters,P,sig,interl(:,Xcshat_int'),deinterl(:,Xcshat_int'));
        
%------Book-keeping for Errors--------
   for u=1:Ka
     num_err=num_err+(~isempty(find(xhat(:,u)-x(:,u),1)));
   end
   %fprintf('num_err=%d,Es/N0=%3.2f dB,Eb/N0=%3.2f\n',num_err,EsN0_vec(sigIter),EbN0_vec(sigIter));    
   berrPrb=num_err/Ka;
   errMat(sigIter,1:4)=[EbN0_vec(sigIter) berrPrb num_err EsN0_vec(sigIter)];
end