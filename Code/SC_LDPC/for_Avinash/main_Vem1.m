% LASSO flag=0 => LASSO w/o any constraints
% LASSO flag=1 => non negative Least Squares lsqnonneg
clear all
T=100;
collision_flag=0;
Es_N0=-16;
Berr=zeros(length(Es_N0),4); %[All T present, one missing, two missing, more than 2 missing]
listSize=T+10;
numSims=10;
load rand_Four_750cross2pow14;

Ac=sqrt(2)*[real(A);imag(A)];
clear A;
J=size(Ac,1);H=size(Ac,2);


for snrIdx=1:length(Es_N0)
    
    fprintf('T=%d,J=%d,H=%d,Es_N0=%3.2fdB,list=%d\n',T,J,H,Es_N0(snrIdx),listSize);
    fprintf('Using nnLS+list size of %d\n',listSize);
    
    sig=sqrt(1/2/(10^(Es_N0(snrIdx)/10))); P=1;
    
    q=0;
    tot_listSize=0;
    while q<numSims
        X_coeffs=randi(H,[1,T]);
        if ~ collision_flag
            if length(unique(X_coeffs))<T
                continue
            end
        end
        
        
        Y=sum(Ac(:,X_coeffs),2)+(sig)*(randn(J,1));
        
        tic
        Xhat=lsqnonneg(Ac,Y);
        toc
        [sortedXhat,sortIndices]=sort(Xhat);
        
        listfinal=sortIndices(end-listSize+1:end);
        
        if (sortedXhat(end-1)-sortedXhat(end-2))>0.5
            listfinal = [listfinal; sortIndices(end-1);sortIndices(end)];
            
        elseif (sortedXhat(end)-sortedXhat(end-1))>0.5
            listfinal = [listfinal; sortIndices(end)];
        end
        
        listfinal=sort(listfinal);
        if ~collision_flag
            missing_idx=0;
            for i=1:length(X_coeffs)
                if ~ismember(X_coeffs(i),listfinal)
                    missing_idx=missing_idx+1;
                end
            end
        end
        
        switch missing_idx
            case 0
                Berr(snrIdx,1)=Berr(snrIdx,1)+1;
            case 1
                Berr(snrIdx,2)=Berr(snrIdx,2)+1;
            case 2
                Berr(snrIdx,3)=Berr(snrIdx,3)+1;
            otherwise
                Berr(snrIdx,4)=Berr(snrIdx,4)+1;
        end
        q=q+1;
        fprintf('TotalSims:%d,e100:%d,e99:%d,e98:%d,eMore:%d\n',q,Berr(snrIdx,1),Berr(snrIdx,2),Berr(snrIdx,3),Berr(snrIdx,4));
    end
end
%%-------------------Results---------------------
%{ 
H=2^14; J=2100; collision=0
Ka=100,Kb=10
e0-denotes ops with 0 indices in error
e1-denotes ops with 1 index in error
 EsN0   sims e0 e1 e2 emore   EbN0
-15.00  10   10              EsN0+10*log10(2100/14)
-17.00  10    9  1           EsN0+10*log10(2100/14)
-18.00  10    0  5  2  3     EsN0+10*log10(2100/14)

H=2^14; J=1500; collision=0
Ka=100,Kb=10
 EsN0   sims e0 e1 e2 emore EbN0
-16.00  10    4  2  3  1    EsN0+10*log10(1500/14)


H=2^14; J=2100; collision=0
Ka=75, Kb=8
 EsN0  sims e0 e1 e2 emore 
-19.00  10  0  0  0  10
-18.00  3   0  1  2  0

%}