function [Vcon_out Ccon1_out Ccon2_out isFix2Cycles]=fix2cycles(Vcon,Ccon1,Ccon2,dv1,dv2)

isFix2Cycles=1;
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

% Fix 2-Cycles Ccon1
for i=1:n
nbhd=Vcon_C(i,1:dv1);
nbhd2=sort(nbhd(nbhd>0),'ascend');
dnbhd=diff(nbhd2); d=sum(dnbhd==0);
while(d>0)
    isFix2Cycles=0;
    rep_Ccon1_ind=nbhd2(1+find(dnbhd==0,1,'first'));
    rand_Ccon1_ind=randi(m1,1,1);
    new_Vcon_ind=Ccon1_V(rand_Ccon1_ind,1);

    Vcon_i_nbhd_ind=find(Vcon_C(i,1:dv1)==rep_Ccon1_ind,1,'first');
    Vcon_new_nbhd_ind=find(Vcon_C(new_Vcon_ind,1:dv1)==rand_Ccon1_ind,1,'first');
    Ccon1_rep_nbhd_ind=find(Ccon1_V(rep_Ccon1_ind,:)==i,1,'first');
    Ccon1_rand_nbhd_ind=1;
    
    Vcon_C(i,Vcon_i_nbhd_ind)=rand_Ccon1_ind;
    Vcon_C(new_Vcon_ind,Vcon_new_nbhd_ind)=rep_Ccon1_ind;
    Ccon1_V(rep_Ccon1_ind,Ccon1_rep_nbhd_ind)=new_Vcon_ind;
    Ccon1_V(rand_Ccon1_ind,Ccon1_rand_nbhd_ind)=i;
    
    s1=Ccon1(rep_Ccon1_ind,Ccon1_rep_nbhd_ind);
    s2=Ccon1(rand_Ccon1_ind,Ccon1_rand_nbhd_ind);
    
    Ccon1_out(rep_Ccon1_ind,Ccon1_rep_nbhd_ind)=s2;
    Ccon1_out(rand_Ccon1_ind,Ccon1_rand_nbhd_ind)=s1;

    nbhd=Vcon_C(i,1:dv1);
    nbhd2=sort(nbhd(nbhd>0),'ascend');
    dnbhd=diff(nbhd2); d=sum(dnbhd==0);
end
end

% Fix 2-Cycles Ccon2
for i=1:n
nbhd=Vcon_C(i,dv1+1:dv);
nbhd2=sort(nbhd(nbhd>0),'ascend');
dnbhd=diff(nbhd2); d=sum(dnbhd==0);
while(d>0)
    isFix2Cycles=0;
    rep_Ccon2_ind=nbhd2(1+find(dnbhd==0,1,'first'));
    rand_Ccon2_ind=randi(m,1,1);
    new_Vcon_ind=Ccon2_V(rand_Ccon2_ind,1);

    Vcon_i_nbhd_ind=dv1+find(Vcon_C(i,dv1+1:dv)==rep_Ccon2_ind,1,'first');
    Vcon_new_nbhd_ind=dv1+find(Vcon_C(new_Vcon_ind,dv1+1:dv)==rand_Ccon2_ind,1,'first');
    Ccon2_rep_nbhd_ind=find(Ccon2_V(rep_Ccon2_ind,:)==i,1,'first');
    Ccon2_rand_nbhd_ind=1;
    
    Vcon_C(i,Vcon_i_nbhd_ind)=rand_Ccon2_ind;
    Vcon_C(new_Vcon_ind,Vcon_new_nbhd_ind)=rep_Ccon2_ind;
    Ccon2_V(rep_Ccon2_ind,Ccon2_rep_nbhd_ind)=new_Vcon_ind;
    Ccon2_V(rand_Ccon2_ind,Ccon2_rand_nbhd_ind)=i;
    
    s1=Ccon2_out(rep_Ccon2_ind,Ccon2_rep_nbhd_ind);
    s2=Ccon2_out(rand_Ccon2_ind,Ccon2_rand_nbhd_ind);
    
    Ccon2_out(rep_Ccon2_ind,Ccon2_rep_nbhd_ind)=s2;
    Ccon2_out(rand_Ccon2_ind,Ccon2_rand_nbhd_ind)=s1;

    nbhd=Vcon_C(i,dv1+1:dv1+dv2);
    nbhd2=sort(nbhd(nbhd>0),'ascend');
    dnbhd=diff(nbhd2); d=sum(dnbhd==0);
end
end

fprintf('Done with fixing 2-cycles...\n');
