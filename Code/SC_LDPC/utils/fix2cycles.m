function [Vcon_out Ccon_out isFix2Cycles]=fix2cycles(Vcon,Ccon)

isFix2Cycles=1;
Vcon_out=Vcon;
Ccon_out=Ccon;

e=max(max(max(Vcon(:)),max(Ccon(:))));
dv=size(Vcon,2); dc=size(Ccon,2);
n=size(Vcon,1); m=size(Ccon,1);
%assert(dv==dv1+dv2);

E=zeros(1,e);
E(Ccon)=repmat((1:m)',1,dc);  E(e)=0;
Vcon_C=E(Vcon);

E=zeros(1,e);
E(Vcon)=repmat((1:n)',1,dv); E(e)=0;
Ccon_V=E(Ccon);

% Fix 2-Cycles Ccon1
for i=1:n
nbhd=Vcon_C(i,1:dv);
nbhd2=sort(nbhd(nbhd>0),'ascend');
dnbhd=diff(nbhd2); d=sum(dnbhd==0);
while(d>0)
    isFix2Cycles=0;
    rep_Ccon_ind=nbhd2(1+find(dnbhd==0,1,'first'));
    rand_Ccon_ind=randi(m,1,1);
    new_Vcon_ind=Ccon_V(rand_Ccon_ind,1);

    Vcon_i_nbhd_ind=find(Vcon_C(i,1:dv)==rep_Ccon_ind,1,'first');
    Vcon_new_nbhd_ind=find(Vcon_C(new_Vcon_ind,1:dv)==rand_Ccon_ind,1,'first');
    Ccon_rep_nbhd_ind=find(Ccon_V(rep_Ccon_ind,:)==i,1,'first');
    Ccon_rand_nbhd_ind=1;
    
    Vcon_C(i,Vcon_i_nbhd_ind)=rand_Ccon_ind;
    Vcon_C(new_Vcon_ind,Vcon_new_nbhd_ind)=rep_Ccon_ind;
    Ccon_V(rep_Ccon_ind,Ccon_rep_nbhd_ind)=new_Vcon_ind;
    Ccon_V(rand_Ccon_ind,Ccon_rand_nbhd_ind)=i;
    
    s1=Ccon(rep_Ccon_ind,Ccon_rep_nbhd_ind);
    s2=Ccon(rand_Ccon_ind,Ccon_rand_nbhd_ind);
    
    Ccon_out(rep_Ccon_ind,Ccon_rep_nbhd_ind)=s2;
    Ccon_out(rand_Ccon_ind,Ccon_rand_nbhd_ind)=s1;
    
    nbhd=Vcon_C(i,1:dv);
    nbhd2=sort(nbhd(nbhd>0),'ascend');
    dnbhd=diff(nbhd2); d=sum(dnbhd==0);
end
end

fprintf('Done with fixing 2-cycles...\n');
