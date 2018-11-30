clear all
addpath('~/Github/avinash/ResearchRep/Code/SC_LDPC/for_Avinash/');
T=125;
maxSNR=3.75;
EsN0_vec=[-14.5 -14];
listSize=T+ceil(T/10);
numSims=15;
load rand_Four_750cross2pow14;

Ac=sqrt(2)*[real(A);imag(A)];
clear A;
J=size(Ac,1);H=size(Ac,2);

fprintf('T=%d,J=%d,H=%d. Using nnLS+list size of %d\n',T,J,H,listSize);

load('chcode_results.mat')
errMat=Hash([num2str(T) ',' num2str(listSize-T) ',0,1']);

idx=find(errMat(:,1)<0,1);
if ~isempty(idx)
    errMat=errMat(1:idx-1,:);
end
Eb_N01=errMat(:,1);

for snrIdx=1:length(EsN0_vec)
    num_err=zeros(length(Eb_N01),1);
    Eb_N0=EsN0_vec(snrIdx)+10*log10(J/log2(H));
    snrEq=reshape(snr_combining(Eb_N01,Eb_N0),[],1);
    
    fprintf('Es_N0=%3.2fdB,Eb_N0=%3.2f\n',EsN0_vec(snrIdx),Eb_N0);
    fprintf('Equiv Eb/N0 is');
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
    while q<numSims 
        X_coeffs=sort(randi(H,[1,T]));
        repeatV=repeat_values(X_coeffs);
        uniq_coeffs=unique(X_coeffs);
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
        
        if (sortedXhat(end-1)-sortedXhat(end-2))>0.5
            listfinal = [listfinal; sortIndices(end-1);sortIndices(end)];
            
        elseif (sortedXhat(end)-sortedXhat(end-1))>0.5
            listfinal = [listfinal; sortIndices(end)];
        end
        
        listfinal=sort(listfinal);
        output_repeatV=repeat_values(listfinal);
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
        if missing_idx<0
            fprintf('WRONG');
        elseif missing_idx>2 || collision_flag<0 || collision_flag>2
            num_err=num_err+T*ones(length(Eb_N01),1);
        else    
            keyvalue=[num2str(T) ',' num2str(listSize-T) ',' num2str(missing_idx) ',' num2str(collision_flag)];
            errMat=Hash(keyvalue);
            num_err=num_err+T*errMat(:,2);
        end
        q=q+1;    
        if mod(q,2)==0
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
%-----------Results--------------
%{
 -------Ka=25--------------
 Eb_N0_eq   biterrPrb Eb_N0_CS
    1.6496    0.7000    4.2996
    2.0261    0.4780    4.2996
    2.4130    0.1220    4.2996
    2.8097    0.0517    4.2996
    3.2155    0.0253    4.2996
    3.6298    0.0199    4.2996
    4.0520    0.0160    4.2996

-------Ka=50--------------
Eb_N0eq berr  snr_CS snr_ch
1.78 7.07e-01 4.80 1.01 
2.15 5.05e-01 4.80 1.51 
2.53 3.60e-01 4.80 2.01 
2.91 6.29e-02 4.80 2.51 
3.31 4.20e-02 4.80 3.01 
3.72 1.27e-02 4.80 3.51 
4.13 8.71e-03 4.80 4.01 
coll stats are 20,0,0,0 
missing_idx stats are 16,1,3,0 

-------Ka=75--------------
coll stats are 15,4,1,0 
missing_idx stats are 14,3,1,2 
    Eb_N0eq berr  snr_CS snr_ch
    1.7840    0.8713    4.7996    1.0103
    2.1496    0.6680    4.7996    1.5103
    2.5261    0.2940    4.7996    2.0103
    2.9130    0.1714    4.7996    2.5103
    3.3097    0.1313    4.7996    3.0103

coll stats are 18,2,0,0 
missing_idx stats are 18,1,1,0 
Eb_N0eq berr  snr_CS snr_ch
2.09 8.43e-01 5.80 1.01 
2.43 6.19e-01 5.80 1.51 
2.78 1.90e-01 5.80 2.01 
3.15 6.82e-02 5.80 2.51 
3.53 2.77e-02 5.80 3.01 
3.91 -7.94e-01 5.80 3.51 
4.31 -7.95e-01 5.80 4.01 

-------Ka=100--------------
coll stats are 13,7,0,0 
missing_idx stats are 11,8,1,0 
Eb_N0eq berr  snr_CS snr_ch
1.78 9.53e-01 4.80 1.01 
2.15 7.47e-01 4.80 1.51 
2.53 2.94e-01 4.80 2.01 
2.91 1.28e-01 4.80 2.51 
3.31 4.41e-02 4.80 3.01 
3.72 2.46e-02 4.80 3.51 
4.13 1.63e-02 4.80 4.01 

coll stats are 16,4,0,0 
missing_idx stats are 17,2,1,0 
Eb_N0eq berr  snr_CS snr_ch
1.93 9.55e-01 5.30 1.01 
2.28 7.50e-01 5.30 1.51 
2.65 2.27e-01 5.30 2.01 
3.03 1.08e-01 5.30 2.51 
3.41 2.83e-02 5.30 3.01 
3.81 1.79e-02 5.30 3.51 
4.22 1.26e-02 5.30 4.01

-------Ka=125--------------

coll stats are 7,5,0,0 
missing_idx stats are 11,1,0,0 
Eb_N0eq berr  snr_CS snr_ch
2.09 9.95e-01 5.80 1.01 
2.43 9.53e-01 5.80 1.51 
2.78 4.78e-01 5.80 2.01 
3.15 2.02e-01 5.80 2.51 
3.53 1.35e-01 5.80 3.01 
3.91 9.42e-02 5.80 3.51 
4.31 9.39e-02 5.80 4.01 

coll stats are 11,3,0,0 
missing_idx stats are 12,2,0,0 
Eb_N0eq berr snr_CS snr_ch
2.26 9.98e-01 6.30 1.01 
2.59 9.56e-01 6.30 1.51 
2.93 4.03e-01 6.30 2.01 
3.28 1.25e-01 6.30 2.51 
3.65 6.16e-02 6.30 3.01 
4.03 1.83e-02 6.30 3.51 
4.41 1.63e-02 6.30 4.01 
%}

