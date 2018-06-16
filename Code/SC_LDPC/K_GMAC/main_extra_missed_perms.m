%{
 As the title says, we consider the scenario where a
 list of Ka-Km correct, Km missed and Kb extra interleavers given as input to decoder.
%}
function errMat=main_extra_missed_perms(Kparams,params,chcode)
global M repeat_flag repeat interl deinterl

N=chcode.N;
M=chcode.M;
L=chcode.L; R=chcode.R;
base_rate=chcode.base_rate;
repeat_flag=chcode.repeat_flag;
repeat=chcode.repeat;

Vcon=chcode.Vcon;
Ccon=chcode.Ccon;
eMax=chcode.eMax;

Ka=Kparams.Ka;
Kb=Kparams.Kb;
Km=Kparams.Km;

collision_flag=params.collision_flag;
if collision_flag
    missing_idx=params.missing_idx;
end
EsN0_vec=params.EsN0_vec;


maxIters=500;
maxSims=20;
max_numerr=max(10,maxSims*Ka*0.01);

lavg=sum(L(1,:).*L(2,:));
ravg=sum(R(1,:).*R(2,:));
assert(int32(10000*ravg*(1-base_rate))==int32(10000*lavg));
%{
fprintf('L(x)=');
for i=1:numel(L(1,:))
   fprintf('%3.2f x^%d+',L(1,i),L(2,i));
end
fprintf('\n R(x)=');
for i=1:numel(R(1,:))
   fprintf('%3.2f x^%d+',R(1,i),R(2,i));
end
fprintf('\n');
%}
n=length(Vcon(:,1)); dv=size(Vcon,2);
m=length(Ccon(:,1)); %Number of check nodes
interl=zeros(N,Ka);

EbN0_vec=EsN0_vec+10*log10(n/(n-m));
sigVec=sqrt(0.5./10.^(EsN0_vec./10));
errMat=-1*ones(length(EbN0_vec),8); 
%-------------Simulations-----------------------------------------
fprintf('Sys.params: N=%d,Ka=%d, (k,n)=(%d,%d) \n',N,Ka,n-m,n);

for i=1:Ka+Kb
    interl(:,i)=reshape(randperm(N),[],1);
    deinterl(interl(:,i),   i)=1:N;
end

if collision_flag>0  %collisionflag=1 implies one double collision
    interl(:,1)=interl(:,round(Ka/2)); 
    deinterl(interl(:,1),1)=1:N;
   if collision_flag>1 %collisionflag=2 implies two double collision
    interl(:,2)=interl(:,round(Ka/2)+1); 
    deinterl(interl(:,2),2)=1:N;
   end
end

chan='Gaussian';

for sigIter=1:length(sigVec)
    sig=sigVec(sigIter);
    fprintf('Es/N0=%3.2f dB,Eb/N0=%3.2f\n',EsN0_vec(sigIter),EbN0_vec(sigIter));

    num_err=0;  sim_cnt=0;
    tic
    while num_err<max_numerr && sim_cnt<maxSims
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
    
    for i=Ka-Km+1:Ka
       interl(:,i)=reshape(randperm(N),[],1);
       deinterl(interl(:,i),i)=1:N;
    end
     %{
   Old logic, where missing index is either one of the collidingindex or
   completely random one. The choice depends on the input missing_idx
    if collision_flag 
        
        for j=1:numel(missing_idx)
          interl(:,missing_idx(j))=reshape(randperm(N),[],1);
          deinterl(interl(:,missing_idx(j)),missing_idx(j))=1:N;
        end
        if Km>numel(missing_idx)
            Km=Km-numel(missing_idx);
            for i=Ka-Km+1:Ka
              interl(:,i)=reshape(randperm(N),[],1);
              deinterl(interl(:,i),i)=1:N;
            end
        end
    else
       for i=Ka-Km+1:Ka
         interl(:,i)=reshape(randperm(N),[],1);
         deinterl(interl(:,i),i)=1:N;
       end
    end
     %}
     indianjasmin v
 ,
    y=sum(x,2)+sig*randn(N,1);
    xhat=BP_KGMAC_extraperms(y,Vcon,Ccon,eMax,maxIters,P,sig,interl,deinterl);
        
%------Book-keeping for Errors--------
fprintf('\n q=%d    :',sim_cnt);
       for u=1:Ka
        err_thissim=numel(find(xhat(:,u)-x(:,u)));
         if err_thissim>0
            num_err=num_err+1;
            fprintf('%d-%d,',u,err_thissim)
         end
       end
       
       
       sim_cnt=sim_cnt+1;
       
       if mod(sim_cnt,2)==0
           fprintf('#err users/total= %d/%d---%3.2f sec\n',num_err,Ka*sim_cnt,toc);
           berrPrb=num_err/sim_cnt/Ka;
        if berrPrb<0.005 && num_err>0   
             break
        end
       end
    end
    
    fprintf('Done. #users_in_err/total= %d/%d---%3.2f sec\n',num_err,Ka*sim_cnt,toc);
    berrPrb=num_err/sim_cnt/Ka;
   if berrPrb<0.005 && num_err>0   
         break
     end
    errMat(sigIter,1:8)=[EbN0_vec(sigIter) berrPrb num_err sim_cnt EsN0_vec(sigIter) Ka Kb Km];
end