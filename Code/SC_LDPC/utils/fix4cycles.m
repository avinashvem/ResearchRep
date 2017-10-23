function [Vcon_out Ccon_out isFix4Cycles]=fix4cycles(Vcon,Ccon)

isFix4Cycles=1;
Vcon_out=Vcon;
Ccon_out=Ccon;

e=max(max(max(Vcon(:)),max(Ccon(:))));
dv=size(Vcon,2); dc=size(Ccon,2);
n=size(Vcon,1); m=size(Ccon,1);
% assert(dv==dv+dv2);

E=zeros(1,e);
E(Ccon)=repmat((1:m)',1,dc);
if numel(find(Ccon==e))>1
    E(e)=0;
end
Vcon_C=E(Vcon);

E=zeros(1,e);
E(Vcon)=repmat((1:n)',1,dv);
if numel(find(Vcon==e))>1
    E(e)=0;
end
Ccon_V=E(Ccon);

for i=1:n
    nbhd=Vcon_C(i,1:dv);
    nbhd=nbhd(nbhd>0);
    i_done=0;
    for j=1:numel(nbhd)
        if(i_done)
            break;
        end
        nbhdj=Ccon_V(nbhd(j),:); nbhdj=nbhdj(nbhdj>0);
        for k=(j+1):numel(nbhd)
            if(i_done)
                break;
            end
            nbhdk=Ccon_V(nbhd(k),:); nbhdk=nbhdk(nbhdk>0);
            com_nbhd=intersect(nbhdj,nbhdk);
            com_nbhd(com_nbhd==i)=[];
            if(numel(com_nbhd)>0)
                i_done=1;
                isFix4Cycles=0;
                k_Ccon_ind=nbhd(k);
                rand_Ccon_ind=randi(m,1,1);
                new_Vcon_ind=Ccon_V(rand_Ccon_ind,1);
                
                Vcon_i_nbhd_k_ind=find(Vcon_C(i,1:dv)==k_Ccon_ind,1,'first');
                Vcon_new_nbhd_rand_ind=find(Vcon_C(new_Vcon_ind,1:dv)==rand_Ccon_ind,1,'first');
                Ccon_k_nbhd_i_ind=find(Ccon_V(k_Ccon_ind,:)==i,1,'first');
                Ccon_rand_nbhd_new_ind=1;
                
                Vcon_C(i,Vcon_i_nbhd_k_ind)=rand_Ccon_ind;
                Vcon_C(new_Vcon_ind,Vcon_new_nbhd_rand_ind)=k_Ccon_ind;
                Ccon_V(k_Ccon_ind,Ccon_k_nbhd_i_ind)=new_Vcon_ind;
                Ccon_V(rand_Ccon_ind,Ccon_rand_nbhd_new_ind)=i;
        
                s1=Ccon(k_Ccon_ind,Ccon_k_nbhd_i_ind);
                s2=Ccon(rand_Ccon_ind,Ccon_rand_nbhd_new_ind);
        
                Ccon_out(k_Ccon_ind,Ccon_k_nbhd_i_ind)=s2;
                try Ccon_out(rand_Ccon_ind,Ccon_rand_nbhd_new_ind)=s1;
                catch
                    break
                end
            end
        end
    end
end

fprintf('Done with fixing 4-cycles...\n');
