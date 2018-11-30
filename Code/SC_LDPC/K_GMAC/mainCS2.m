%{
Contrary to mainCS, here for each simulation we call the channel coding
part after CS decoder is finished. Output from CS decoder is input to ch. coding.
%}
clear all
addpath('~/Github/avinash/ResearchRep/Code/SC_LDPC/for_Avinash/');

%% ---Channel code initialization-----------
chcode.N=28500;
chcode.M=144; %(144,84) code
chcode.L=[0.5 0.5;2 3];chcode.R=[1;6];
chcode.base_rate=7/12;
chcode.repeat_flag=0;
chcode.repeat=1;
girth=0;
%--Michael's PEG construction------------------
while girth<6
[chcode.Vcon,chcode.Ccon,chcode.eMax,girth]=LDPC_PEG_cnstr(chcode.M,round(chcode.M*(1-chcode.base_rate)),chcode.L,chcode.R,chcode.repeat_flag);
end

%% -------------
T=175;
EsN0_vec=[-12.75 -12.25 -11.75];
Es_N01=[2.75 3.25 4];
maxSNR=6;
numSims=15;

Eb_N01=Es_N01 + 10*log10(1/chcode.base_rate);
listSize=T+ceil(T/10);

load rand_Four_750cross2pow14;
Ac=sqrt(2)*[real(A);imag(A)];
clear A;
J=size(Ac,1);H=size(Ac,2);

fprintf('T=%d,J=%d,H=%d. Using nnLS+list size of %d\n',T,J,H,listSize);

for snrIdx=1:length(EsN0_vec)
    num_err=zeros(length(Eb_N01),1);
    Eb_N0=EsN0_vec(snrIdx)+10*log10(J/log2(H));
    
    snrEq=snr_combining(Eb_N01,Eb_N0);
    
    fprintf('Es_N0=%3.2fdB,Eb_N0=%3.2f\n',EsN0_vec(snrIdx),Eb_N0);
    fprintf('Equiv Eb/N0 is ');
    for zz=1:numel(snrEq)
        fprintf('%3.2f,',snrEq(zz));
    end
    fprintf('\n');
    
    sig=sqrt(1/2/(10^(EsN0_vec(snrIdx)/10))); P=1;
    
    q=0;
    tot_listSize=0;
    Mat_coll=[0 0 0 0]; %fraction of collision_flag being 0 1 2 -1
    Mat_Km=[0 0 0 0]; %fraction of missing_idx being 0 1 2 more
    
    tic
    while q<numSims && min(num_err)<1.5*(0.05*numSims*T)
        X_coeffs=sort(randi(H,[1,T]));
        
        repeatV=repeat_values(X_coeffs);
        uniq_coeffs=unique(X_coeffs);
        
   % Just to keep track of collision statistics
        switch length(uniq_coeffs)
            case T
                collision_flag=0;
            case T-1
                collision_flag=1;
            case T-2
                sort_Xcoeffs=sort(X_coeffs);
                diff=[sort_Xcoeffs sort_Xcoeffs(end)+1]-[0 sort_Xcoeffs];
                zero_idx=find(diff==0);
                if numel(zero_idx)~=2
                    fprintf('debug')
                elseif zero_idx(2)==zero_idx(1)+1 %Triple collision
                    collision_flag=-1;
                else
                    collision_flag=2; %Two double collision
                end
            otherwise 
                collision_flag=-1;
        end
        
        Y=sum(Ac(:,X_coeffs),2)+(sig)*(randn(J,1));
        Xhat=lsqnonneg(Ac,Y);
        
        [sortedXhat,sortIndices]=sort(Xhat);
        listfinal=sortIndices(end-listSize+1:end);
% Accounting for repeated indices in decoder     
        if (sortedXhat(end-1)-sortedXhat(end-2))>0.5
            listfinal = [listfinal; sortIndices(end-1);sortIndices(end)];            
        elseif (sortedXhat(end)-sortedXhat(end-1))>0.5
            listfinal = [listfinal; sortIndices(end)];
        end
        
        listfinal=sort(listfinal);
        output_repeatV=repeat_values(listfinal);
% The whole missing_idx biz is unnecessary for this
        missing_idx=numel(setdiff(X_coeffs,listfinal));  
        if ~isempty(repeatV)
            for i=1:numel(repeatV(:,1))
                if ~isempty(output_repeatV)
                    if ismember(repeatV(i,1),output_repeatV(:,1))
                        Idx=find(output_repeatV(:,1)==repeatV(i,1));
                        missing_idx=missing_idx+max(repeatV(Idx,2)-output_repeatV(Idx,2),0);
                    end
                else
                    missing_idx=missing_idx+sum(repeatV(:,2)-1);
                end                
            end
        end
% Keeping track of missing_idx and collision statistics        
        if missing_idx>2
            Mat_Km(4)=Mat_Km(4)+1;
        else
            Mat_Km(missing_idx+1)=Mat_Km(missing_idx+1)+1;
        end
        if collision_flag==-1
            Mat_coll(4)=Mat_coll(4)+1;
        else
            Mat_coll(collision_flag+1)=Mat_coll(collision_flag+1)+1;
        end
            
%  fprintf('Km=%d,Coll=%d\n',missing_idx,collision_flag);
        params.Ka=T;
        params.EsN0_vec=Es_N01;
        errMat=main_extra_missed_perms1(chcode,params,X_coeffs,listfinal);
        num_err=num_err+reshape(errMat(:,3),length(Es_N01),1);
        q=q+1;    
        if mod(q,1)==0
            fprintf('q=%d,num_err is',q);
            for zz=1:numel(num_err)
                fprintf('%d,',num_err(zz));
            end
            fprintf(',, In %4.3f secs.\n',toc);
            fprintf('coll stats are %d,%d,%d,%d \n',Mat_coll);
            fprintf('missing_idx stats are %d,%d,%d,%d \n',Mat_Km);
        end
        
         upperlim_SNR=find(snrEq(:,1)>maxSNR,1);
        if ~isempty(upperlim_SNR) 
            if num_err(upperlim_SNR) >1.5*(0.05*numSims*T)
                fprintf('For eq SNR of %3.2f, errs is %d and errPrb for 20 sims=%4.3e\n',snrEq(upperlim_SNR),num_err(upperlim_SNR),num_err(upperlim_SNR)/q/T);
                break
            end
        end
    end
    snrEq=reshape(snr_combining(Eb_N01,Eb_N0),[],1);
    snrEq(:,2)=num_err./(q*T);
    snrEq(:,3)=Eb_N0*ones(length(Eb_N01),1);
    snrEq(:,4)=Eb_N01;
    
    load('final_results.mat');
    keyvalue=num2str(T);
    L=0;
    if HashCS.isKey(keyvalue)
        res=HashCS(keyvalue);
        L=length(res(:,1,1));
    end
    res(L+1).data=snrEq;
    res(L+1).cs_coll=Mat_coll;%0 1 2 -1
    res(L+1).missing=Mat_Km; %0 1 2 more
    HashCS(keyvalue)=res;
    
    save('final_results.mat','HashCS');
    
for i=1:numel(snrEq(:,1))
        fprintf('%3.2f %3.2e %3.2f %3.2f \n',snrEq(i,1),snrEq(i,2),snrEq(i,3),Eb_N01(i));
end
end
%-----------Results via second method
%{
 ----------Ka=25---------
  Eb_N0_eq   biterrPrb Eb_N0_CS
    2.7155    0.1067    3.7996
    3.3400    0.0467    3.7996
    3.9814    0.0400    3.7996
    4.6379    0.0400    3.7996
    5.3077    0.0400    3.7996
Eb_N0_eq   biterrPrb Eb_N0_CS
    2.8097    0.0600    4.2996
    3.4217    0.0360    4.2996
    4.0520    0.0240    4.2996
    4.6987    0.0240    4.2996
    5.3598    0.0240    4.2996
  Eb_N0_eq    biterr   Eb_N0_CS  Eb_N0_ch
    2.7183    0.1000    4.7996    2.2603
    3.1102    0.0533    4.7996    2.7603
    3.5115    0.0200    4.7996    3.2603

----------------Ka=150-----------------
coll stats are 7,2,5,0 
missing_idx stats are 11,1,0,2 
Eb_N0_eq biterrPrb Eb_N0_CS
3.49 2.87e-01 7.00 2.51 
4.02 6.36e-02 7.00 3.26 
4.57 3.60e-02 7.00 4.01 
5.15 3.87e-02 7.00 4.76 
5.75 2.67e-02 7.00 5.51 

----------------Ka=175-----------------

chcode.L=[0.5 0.5;2 3];chcode.R=[1;6]; (144,84) code
coll stats are 7,4,3,1 
missing_idx stats are 12,1,0,2 
5.53 1.32e-01 7.55 5.09 
5.92 8.95e-02 7.55 5.59 
6.53 5.22e-02 7.55 6.34 

----------------Ka=200-----------------

coll stats are 0,0,2,0 
missing_idx stats are 2,0,0,0 
4.59 7.92e-01 8.30 3.51 
5.11 8.30e-01 8.30 4.26 
5.65 7.72e-01 8.30 5.01 
6.22 8.85e-01 8.30 5.76 
6.81 8.90e-01 8.30 6.51 
7.42 8.90e-01 8.30 7.26 


Incomplete;q=2
coll stats are 0,0,2,0 
missing_idx stats are 0,0,0,2 
4.76 6.67e-01 8.80 3.51 
5.26 6.75e-01 8.80 4.26 
5.78 6.52e-01 8.80 5.01 
6.34 7.80e-01 8.80 5.76 
6.91 6.93e-01 8.80 6.51 
7.51 4.70e-01 8.80 7.26 

%}
