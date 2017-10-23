function [index_vn,index_cn,girth]=ldpc_peg_const(n,m,d)
%n: block length of ldpc
%m: # cns
%d: degree of vn
%index_vn: id of cns that each vn connects to
%index_cn: id of vns that each cn conncects to
%------------------------------------------------------------------
deg_vn=zeros(1,n);   % degree of each VN
deg_cn=zeros(1,m);   % degree of each CN
cn_id=zeros(1,m); 
vn_id=zeros(1,n);

firstgirth=1;
girth=0;
ncycle=0;
gnew=0;
for i=1:n
   
    for k=1:d(i)
        if k==1
            noCtree=1:m;
            [num,dc_min_list]=findmin(noCtree,deg_cn);  %need to mod: id is not unique      
            dc_min_index=unidrnd(num);
            dc_min_id=dc_min_list(dc_min_index);
            deg_vn(i)=deg_vn(i)+1;         
            deg_cn(dc_min_id)=deg_cn(dc_min_id)+1;
            index_cn(i,deg_vn(i))=dc_min_id;
            index_vn(dc_min_id,deg_cn(dc_min_id))=i;
            
        else
            l=1;   % #level
            Ctree=zeros(1,m);    %# list of CN in the current tree
            Vtree=zeros(1,n);   %# list of CN in the current tree
            cn_id=zeros(1,m);
            vn_id=zeros(1,n);
            cn_id(1)=deg_vn(i);  % # of CN in the 1st level
            vn_id(1)=1;            %# of VN in the 1st level(root)
            nctree=deg_vn(i);    % # of CN in the current tree
            nvtree=1;               %# of CN in the current tree
            CN(1,1:cn_id(1))=index_cn(i,1:cn_id(1)); %list of CN in the 1st level
            VN(1,1)=i;
            Ctree(1:nctree)=CN(1,1:cn_id(1));
            Vtree(1)=i;
            end_tree_flag1=0;
            end_tree_flag2=0;
            while end_tree_flag1==0 && end_tree_flag2==0              
                l=l+1;              
                for j=1:cn_id(l-1)  %j-index_vn of CN in the l_th level; cn_id(l)-# of CN in the lth level
                    for t=1:deg_cn(CN(l-1,j)) %t-index of VN in the l+1_th level;
                      flag=0;                  %  CN(l,j)-jth CN in the lth level;
                                               %  deg_cn(CN(l,j))-deg of jth
                                               %   CN in the lth level
                      for p=1:nvtree          % index(j,t)-t-th VN connected to CN j
                        if index_vn(CN(l-1,j),t)==Vtree(p)
                            flag=1;
                            break
                        end
                       end
                       if flag==0
                           vn_id(l)=vn_id(l)+1;
                           VN(l,vn_id(l))=index_vn(CN(l-1,j),t);
                           nvtree=nvtree+1;
                           Vtree(nvtree)=index_vn(CN(l-1,j),t);
                       end
                    end
                end
                if vn_id(l)==0 
                    end_tree_flag1=1;
                    break
                end
                for j=1:vn_id(l)
                    for t=1:deg_vn(VN(l,j))
                        flag=0;
                        for p=1:nctree
                            if index_cn(VN(l,j),t)==Ctree(p)
                                flag=1;
                                break
                            end
                        end
                        if flag==0
                           cn_id(l)=cn_id(l)+1;
                           CN(l,cn_id(l))=index_cn(VN(l,j),t);
                           nctree=nctree+1;
                           Ctree(nctree)=index_cn(VN(l,j),t);
                       end
                    end
                end
           
                if nctree==m 
                    end_tree_flag2=1;
                    break
                end                
            end  % end while
            
            %-------------- add an edge---------------------
            if end_tree_flag1==1               
                noCtree=1:m;
                for j=1:nctree
                    noCtree(Ctree(j))=0;
                end
                [num,dc_min_list]=findmin(noCtree,deg_cn);
%                 dc_min_index=unidrnd(num);
                
            elseif end_tree_flag2==1
                [num,dc_min_list]=findmin(CN(l,1:cn_id(l)),deg_cn);
                ncycle=ncycle+1;
                girth_new(ncycle)=2*l;
                gnew=2*l;
                if firstgirth
                    firstgirth=0;
                    girth=gnew;
                end
                if gnew<girth
                    girth=gnew;
                end
              
            end
            dc_min_index=unidrnd(num);
            dc_min_id=dc_min_list(dc_min_index);
            deg_vn(i)=deg_vn(i)+1;
            deg_cn(dc_min_id)=deg_cn(dc_min_id)+1;
            index_cn(i,deg_vn(i))=dc_min_id;
            index_vn(dc_min_id,deg_cn(dc_min_id))=i;
                
        end
    end
    
end
  
% for i=1:n
%     for j=1:3
%         H(index_cn(i,j),i)=1;
%     end
% end
        
            