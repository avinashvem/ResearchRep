clear; clc;
% For the LDPC part
dv1=3; % Variable-node degree
dc1=6; % Check-node degree

% For the LDGM Part
dv2=6; % Variable-node degree
dc2=3; % Check-node degree
dv=dv1+dv2;

c=200;
ms=dv2*(c); % Block Length / # of check-nodes in LDGM
ns=dc2*(c); % # of variable-nodes
m1s=ns*(dv1/dc1); % # of check-nodes in LDPC
e1s=ns*dv1; e2s=ns*dv2; es=ns*dv;

L=15; % Chain length of the SC system
w=3; % coupling parameter
Lw=L+w-1; % code-bits system length
ldpc_start=w; ldpc_end=L+w-1; 
ldpc_len=ldpc_end-ldpc_start+1; % ldpc-check system length
m=ms*L;
n=ns*Lw;
m1=m1s*ldpc_len;

e1=Lw*e1s; e2=Lw*e2s; e=Lw*es;

E=zeros(1,e+1); % This is the main object in the sim

Vcon=zeros(n,dv);
Vcon(:,1:dv1)=reshape(1:e1,dv1,[])';
Vcon(:,dv1+1:dv)=e1+reshape(1:e2,dv2,[])';

P1=cell(Lw); P2=cell(Lw);

Ccon1=zeros(m1,dc1); Ccon2=zeros(m,dc2);

for i=1:Lw
    P1{i}=e1s*(i-1)+reshape(randperm(e1s),w,[]);
    P2{i}=e1+e2s*(i-1)+reshape(randperm(e2s),w,[]);
end

for i=ldpc_start:ldpc_end
    Etmp1=zeros(1,e1s);
    for j=1:w
       Ptmp1=P1{i-j+1};
       Etmp1((j-1)*e1s/w+1:j*e1s/w)=Ptmp1(w+1-j,:);
    end
    Etmp1=Etmp1(randperm(e1s));
    Ccon1((i-ldpc_start)*m1s+1:(i-ldpc_start+1)*m1s,:)=reshape(Etmp1,[],dc1);
end

for i=1:L
    Etmp2=zeros(1,e2s);
    for j=1:w
       Ptmp2=P2{i+j-1};
       Etmp2((j-1)*e2s/w+1:j*e2s/w)=Ptmp2(w+1-j,:);
    end
    Etmp2=Etmp2(randperm(e2s));
    Ccon2((i-1)*ms+1:i*ms,:)=reshape(Etmp2,[],dc2);
end

%%%% Delete unused Vcon sockets
E_code=zeros(1,e+1);
E_code(Ccon1)=ones(size(Ccon1)); E_code(Ccon2)=ones(size(Ccon2));
Vcon(E_code(Vcon)==0)=e+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

isFix4Cycles=0;
while(~isFix4Cycles)
    isFix2Cycles=0;
    while(~isFix2Cycles)
        [Vcon Ccon1 Ccon2 isFix2Cycles]=fix2cycles(Vcon,Ccon1,Ccon2,dv1,dv2);
    end
    [Vcon Ccon1 Ccon2 isFix4Cycles]=fix4cycles(Vcon,Ccon1,Ccon2,dv1,dv2);
end

% Changes such that Vcon and Ccon can be input to Irregular MP Decoder
e_new=e+m; %m extra edges between LDGM check nodes and bit nodes from chan
Ccon2(:,dc2+1)=e+(1:m)';

if dc2+1<dc1
    Ccon2(:,dc2+2:dc1)=e_new+1; dc=dc1;
elseif dc2+1>dc1 
    Ccon1(:,dc1+1:dc2+1)=e_new+1; dc=dc2+1;
end

Ccon=[Ccon1;Ccon2];        Ccon1=Ccon1(:,1:dc1); Ccon2=Ccon2(:,1:dc2);
Vcon(Vcon==e+1)=e_new+1;
Vcon_extra=(e_new+1)*ones(m,dv);     Vcon_extra(:,1)=e+(1:m)';
Vcon=[Vcon; Vcon_extra];

fid = fopen('Vcon_CC_Coupled_3_6_6_3_200','w');
fprintf(fid,'Vcon Matrix\n');
fprintf( fid,'%d %d %d %d %d %d %d %d %d \n', Vcon');
fclose(fid);

fid = fopen('Ccon_CC_Coupled_3_6_6_3_200','w');
fprintf(fid,'Ccon Matrix\n');
fprintf(fid,'%d %d %d %d %d %d \n', Ccon');
fclose(fid);

fid = fopen('Simulation_Parameters','w');
fprintf(fid,'%d %d %d %d %d %d %d %d %d %d %d %d\n',dv1,dc1,dv2,dc2,m, m1, n, n+m, m+m1, e_new,L,w);



