clear all
close all
clc
tic

k=200; % number of information bits
n=400; % number of codeword bits
l=n-k; % number of rows in H
col_weight=3;
row_weight=6; % row weight is a constant always.

% assign weight to each variable node
var_deg_vec=col_weight*ones(1,n); % vector which contains degree of each variable node in H. 

% assign weight to each check node
check_deg_vec=row_weight*ones(1,l); % vector which contains degree of each check node in H.   

H=zeros(l,n); % initialize parity check matrix
overall_girth=inf;
lmax=100; % This is the maximum allowable girth, play around with this to see what the speed is 

%================= no need to modify anything below this line=============
% PEG algorithm begins
for jj=1:n
    
    cntr=1;
    progress=jj;
    
    if cntr==1
        
        check_cand=1:l;%possible check candidates
        
        check_cand_deg=sum(H(check_cand,:),2);%degree of check candidates
        
        temp=find(check_cand_deg==min(check_cand_deg));
        
        ptr=ceil(rand*length(temp));
        
        H(check_cand(temp(ptr)),jj)=1;
        
        cntr=cntr+1;
        
        %H
        
    end
    
    while cntr > 1 && cntr<=var_deg_vec(jj)
        
        girth=2;
        
        nbrhd_check=zeros(l,1);%neighborhood in check nodes
        
        curr_nbrhd_check=zeros(l,1);
        
        curr_nbrhd_var=zeros(1,n);%current neighborhood in variable nodes
        
        temp_check=find(H(:,jj)==1);
        
        nbrhd_check(temp_check)=1;
        
        curr_nbrhd_check(temp_check)=1;
        
        tempp_vec=setdiff((1:l),temp_check);
        
        nbrhd_check_prev_comp=tempp_vec;
        
        temp_deg_vec=sum(H(nbrhd_check_prev_comp,:),2)';
        
        nbrhd_check_prev_comp=nbrhd_check_prev_comp(temp_deg_vec<check_deg_vec(nbrhd_check_prev_comp));
        
        l_n_check_prev=sum(nbrhd_check);
        
        indctr=1;
        
        while indctr==1
        
        curr_nbrhd_var=zeros(1,n);
        
        temppp=find(curr_nbrhd_check==1);
        
        for mm=1:sum(curr_nbrhd_check)
            
            curr_nbrhd_var(find(H(temppp(mm),:)==1))=1;
            
        end
        
        curr_nbrhd_check=zeros(l,1);
        
        temppp=find(curr_nbrhd_var==1);
        
        for mm=1:sum(curr_nbrhd_var)
            
            temp=find(H(:,temppp(mm))==1);
            
            curr_nbrhd_check(temp)=1;
            
            nbrhd_check(temp)=1;
            
        end
        
        tempp_vec=setdiff((1:l),find(nbrhd_check==1));
        
        nbrhd_check_curr_comp=tempp_vec;
        
        temp_deg_vec=sum(H(nbrhd_check_curr_comp,:),2)';
        
        nbrhd_check_curr_comp=nbrhd_check_curr_comp(temp_deg_vec<check_deg_vec(nbrhd_check_curr_comp));
        
        l_n_check=sum(nbrhd_check);
        
        if (l_n_check==l_n_check_prev)||(isempty(nbrhd_check_curr_comp)==1)||girth>lmax
            
            %jj
            %girth=girth+2;
            
            if ((isempty(nbrhd_check_curr_comp)==1))
                
                girth=girth+2;
                
            else
                
                girth=inf;
                
            end
            
            %girth
            
            %girth=min(girth,girth_prev);
            
            %girth_prev=girth;
            
            check_cand=nbrhd_check_prev_comp;
            
            check_cand_deg=sum(H(check_cand,:),2);%degree of check candidates
            
            temp=find(check_cand_deg==min(check_cand_deg));
            
            ptr=ceil(rand*length(temp));
            
            H(check_cand(temp(ptr)),jj)=1;
            
            %H
            
            cntr=cntr+1;
            
             indctr=0;
            
        else
            
            l_n_check_prev=l_n_check;
            
            nbrhd_check_prev_comp=nbrhd_check_curr_comp;
            
            girth=girth+2;
            
        end
        
        end
        

    end
    
    overall_girth=min(overall_girth,girth);
    
end

girth=overall_girth
temp=5;

