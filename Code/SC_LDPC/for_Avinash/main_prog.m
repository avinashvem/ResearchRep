% This program is for T=5 i.e T-sparse Compressed sensing problem where the vector is of 0's and 1's.
% LASSO flag=0 => LASSO w/o any constraints
% LASSO flag=1 => non negative Least Squares lsqnonneg

clear all
T=100;
Eb_N0=[-17];
Berr=0;
listSize=T+50;
numFalsePositive=0;
numSims=1;
collisionFlag=1; % 0 for no collision and 1 for collision
%load rand_Four_1050cross2pow14;
load rand_Four_750cross2pow12;

Ac=sqrt(2)*[real(A);imag(A)];
clear A;
J=size(Ac,1);H=size(Ac,2);
permut = randperm(H);
Acpermuted = Ac(:,permut);
%keyboard

csOutListg100=0; 
numErrListg100=0;

csOutListe100=0;
numErrListe100=0;

csOutListe99=0;
numErrListe99=0;

csOutListe98=0;
numErrListe98=0;

numColl=0;

for snrIdx=1:length(Eb_N0)
    
    fprintf('J=%d,H=%d,Es_N0=%3.2fdB,list=%d\n',J,H,Eb_N0(snrIdx),listSize);
    fprintf('Using nnLS+list size of %d\n',listSize);
    
    sig=sqrt(1/2/(10^(Eb_N0(snrIdx)/10))); P=1;
    
    q=0;
    tot_listSize=0;
    %keyboard;
    while q<numSims
        X_coeffs=randi(H,[1,T]);
        if collisionFlag==0
            if length(unique(X_coeffs))~=T
                continue;
            end
        else
            if length(unique(X_coeffs))~=T-1
                continue;
            end
        end
        
        Y=sum(Ac(:,X_coeffs),2)+(sig)*(randn(J,1));
      %  Yprime=sum(Acpermuted(:,X_coeffs),2)+(sig)*(randn(J,1));
        
        tic
        Xhat=lsqnonneg(Ac,Y);
       % Xhatprime=lsqnonneg(Acpermuted,Yprime);
        fprintf('%f secs\n',toc)
        [sortedXhat,sortIndices]=sort(Xhat);
        sortedXhat(end-4:end)'
        %[sortedXhatprime,sortIndicesprime]=sort(Xhatprime);
        %{
        list=sortIndices(end-listSize+1:end);
        listprime=sortIndicesprime(end-listSize+1:end);
        listfinal=intersect(list,listprime);
        if ((sortedXhat(end)+sortedXhatprime(end))/2)-((sortedXhat(end-1)+sortedXhatprime(end-1))/2)>0.4
            numColl=numColl+1
            listfinal = [listfinal; sortIndices(end)];
        end
        L = length(listfinal);
        if length(listfinal)>T
            L=101;
        end
        
        X_C=sort(X_coeffs)';
        LF=sort(listfinal);
        switch L
            case 101
                csOutListg100=csOutListg100+1;
                for i=1:length(X_coeffs)
                    if ~ismember(X_coeffs(i),listfinal)
                        numErrListg100=numErrListg100+1;
                    end
                end
            case 100
                csOutListe100=csOutListe100+1;
                for i=1:length(X_coeffs)
                    if ~ismember(X_coeffs(i),listfinal)
                        numErrListe100=numErrListe100+1;
                    end
                end
            case 99
                csOutListe99=csOutListe99+1;
                for i=1:length(X_coeffs)
                    if ~ismember(X_coeffs(i),listfinal)
                        numErrListe99=numErrListe99+1;
                    end
                end
            case 98
                csOutListe98=csOutListe98+1;
                for i=1:length(X_coeffs)
                    if ~ismember(X_coeffs(i),listfinal)
                        numErrListe98=numErrListe98+1;
                    end
                end
        end
        q=q+1;
        fprintf('g100:%d,%d,e100:%d,%d,e99:%d,%d,e98:%d,%d',csOutListg100,numErrListg100,csOutListe100,numErrListe100,csOutListe99,numErrListe99,csOutListe98,numErrListe98);
        %keyboard
        %}
    end
    
    
end



