function [Vcon_out Ccon1_out Ccon2_out isFix4Cycles]=fix4cycles(Vcon,Ccon1,Ccon2,dv1,dv2)

isFix4Cycles=1;
Vcon_out=Vcon;
Ccon1_out=Ccon1;
Ccon2_out=Ccon2;

e=max(max(max(Vcon(:)),max(Ccon1(:))),max(Ccon2(:)));
dv=size(Vcon,2); dc1=size(Ccon1,2); dc2=size(Ccon2,2);
n=size(Vcon,1); m1=size(Ccon1,1); m=size(Ccon2,1);
assert(dv==dv1+dv2);

E=zeros(1,e);
E(Ccon1)=repmat((1:m1)',1,dc1); E(Ccon2)=repmat((1:m)',1,dc2);
Vcon_C=E(Vcon);

E=zeros(1,e);
E(Vcon)=repmat((1:n)',1,dv);
Ccon1_V=E(Ccon1); Ccon2_V=E(Ccon2);

for i=1:n
    nbhd=Vcon_C(i,1:dv1);
    nbhd=nbhd(nbhd>0);
    i_done=0;
    for j=1:numel(nbhd)
        if(i_done)
            break;
        end
        nbhdj=Ccon1_V(nbhd(j),:); nbhdj=nbhdj(nbhdj>0);
        for k=(j+1):numel(nbhd)
            if(i_done)
                break;
            end
            nbhdk=Ccon1_V(nbhd(k),:); nbhdk=nbhdk(nbhdk>0);
            com_nbhd=intersect(nbhdj,nbhdk);
            com_nbhd(com_nbhd==i)=[];
            if(numel(com_nbhd)>0)
                i_done=1;
                isFix4Cycles=0;
                k_Ccon1_ind=nbhd(k);
                rand_Ccon1_ind=randi(m1,1,1);
                new_Vcon_ind=Ccon1_V(rand_Ccon1_ind,1);
                
                Vcon_i_nbhd_k_ind=find(Vcon_C(i,1:dv1)==k_Ccon1_ind,1,'first');
                Vcon_new_nbhd_rand_ind=find(Vcon_C(new_Vcon_ind,1:dv1)==rand_Ccon1_ind,1,'first');
                Ccon1_k_nbhd_i_ind=find(Ccon1_V(k_Ccon1_ind,:)==i,1,'first');
                Ccon1_rand_nbhd_new_ind=1;
                
                Vcon_C(i,Vcon_i_nbhd_k_ind)=rand_Ccon1_ind;
                Vcon_C(new_Vcon_ind,Vcon_new_nbhd_rand_ind)=k_Ccon1_ind;
                Ccon1_V(k_Ccon1_ind,Ccon1_k_nbhd_i_ind)=new_Vcon_ind;
                Ccon1_V(rand_Ccon1_ind,Ccon1_rand_nbhd_new_ind)=i;
        
                s1=Ccon1(k_Ccon1_ind,Ccon1_k_nbhd_i_ind);
                s2=Ccon1(rand_Ccon1_ind,Ccon1_rand_nbhd_new_ind);
        
                Ccon1_out(k_Ccon1_ind,Ccon1_k_nbhd_i_ind)=s2;
                Ccon1_out(rand_Ccon1_ind,Ccon1_rand_nbhd_new_ind)=s1;
            end
        end
    end
end

for i=1:n
nbhd=Vcon_C(i,dv1+1:dv);
nbhd=nbhd(nbhd>0);
i_done=0;
for j=1:numel(nbhd)
if(i_done)
    break;
end
nbhdj=Ccon2_V(nbhd(j),:); nbhdj=nbhdj(nbhdj>0);
    for k=(j+1):numel(nbhd)
        if(i_done)
            break;
        end
        nbhdk=Ccon2_V(nbhd(k),:); nbhdk=nbhdk(nbhdk>0);
        com_nbhd=intersect(nbhdj,nbhdk);
        com_nbhd(com_nbhd==i)=[];
        if(numel(com_nbhd)>0)
            i_done=1;
            isFix4Cycles=0;
            k_Ccon2_ind=nbhd(k);
            rand_Ccon2_ind=randi(m,1,1);
            new_Vcon_ind=Ccon2_V(rand_Ccon2_ind,1);
            
            Vcon_i_nbhd_k_ind=dv1+find(Vcon_C(i,1+dv1:dv)==k_Ccon2_ind,1,'first');
            Vcon_new_nbhd_rand_ind=dv1+find(Vcon_C(new_Vcon_ind,dv1+1:dv)==rand_Ccon2_ind,1,'first');
            Ccon2_k_nbhd_i_ind=find(Ccon2_V(k_Ccon2_ind,:)==i,1,'first');
            Ccon2_rand_nbhd_new_ind=1;
            
            Vcon_C(i,Vcon_i_nbhd_k_ind)=rand_Ccon2_ind;
            Vcon_C(new_Vcon_ind,Vcon_new_nbhd_rand_ind)=k_Ccon2_ind;
            Ccon2_V(k_Ccon2_ind,Ccon2_k_nbhd_i_ind)=new_Vcon_ind;
            Ccon2_V(rand_Ccon2_ind,Ccon2_rand_nbhd_new_ind)=i;
    
            s1=Ccon2(k_Ccon2_ind,Ccon2_k_nbhd_i_ind);
            s2=Ccon2(rand_Ccon2_ind,Ccon2_rand_nbhd_new_ind);
    
            Ccon2_out(k_Ccon2_ind,Ccon2_k_nbhd_i_ind)=s2;
            Ccon2_out(rand_Ccon2_ind,Ccon2_rand_nbhd_new_ind)=s1;
        end
    end
end
end

fprintf('Done with fixing 4-cycles...\n');
