clear all
addpath('../PEG/')
addpath('~/Github/avinash/ResearchRep/Code/SC_LDPC/for_Avinash/');
chcode.N=28500;
chcode.M=168;
chcode.kmsg=84;
chcode.tBCH=2;
chcode.L=[0.5 0.5;2 3];chcode.R=[1;6];
chcode.base_rate=7/12;
chcode.repeat_flag=0;
chcode.repeat=1;
girth=0;
%--Michael's PEG construction------------------
while girth<6
[chcode.Vcon,chcode.Ccon,chcode.eMax,girth]=LDPC_PEG_cnstr(chcode.M,round(chcode.M*(1-chcode.base_rate)),chcode.L,chcode.R,chcode.repeat_flag);
end
fprintf('[L(x) and R(X) \n');
[chcode.L chcode.R]

Ka_vec=175;
Km_vec=[0 1 2];
collision_vec=[2 1 0];

params.EsN0_vec=[3 3.75 4.5];
params.missing_idx=[];

op=zeros(length(Ka_vec),length(Km_vec),length(collision_vec),length(params.EsN0_vec),8);
for idxK=1:length(Ka_vec)
    Kparams.Ka=Ka_vec(idxK);
    Kparams.Kb=10;%ceil(Kparams.Ka/10);
    
    for idxc=1:length(collision_vec)
        params.collision_flag=collision_vec(idxc);
        
        for idxm=1:length(Km_vec)
            Kparams.Km=Km_vec(idxm);
            load('chcode_results.mat')
            keyvalue=[num2str(Kparams.Ka) ',' num2str(Kparams.Kb) ',' num2str(Kparams.Km) ',' num2str(params.collision_flag)];
%             if Hash.isKey(keyvalue)
%                 fprintf('Already found\n');
%                 continue
%             else
                fprintf('Parameters: Ka=%d,Kb=%d,Km=%d,collision flag=%d\n',Kparams.Ka,Kparams.Kb,Kparams.Km,params.collision_flag)
                errMat=main_extra_missed_perms2(Kparams,params,chcode);
                op(idxK,idxc,idxm,:,:)=errMat;
                
                load('chcode_results.mat')
                Hash(keyvalue)=errMat;
           %    assert(aN==chcode.N);
           %    assert(aM==chcode.M);
           %    save('chcode_results.mat','Hash','aN','aM')
                 
                for i=1:numel(errMat(:,1))
                    fprintf('%3.2f %3.2e %d %d %3.2f \n',errMat(i,1),errMat(i,2),errMat(i,3),errMat(i,4),errMat(i,5));
                end
            %end
          end
      end
for i=1:numel(errMat(:,1))
    fprintf('%3.2f %3.2e %d %d %3.2f \n',errMat(i,1),errMat(i,2),errMat(i,3),errMat(i,4),errMat(i,5));
end
end

%% Results
% Ka=175, Kb=10
% L(x)=x^2/2+x^3/2 ;R(x)=x^6. Rate=7/12=0.5833
% Single code. No outer code. (84,144) base code; 
%{
Km=0; Collision flag=2
Eb_N0   berr   err sims Es_N0
4.59 1.00e-01 35 2 2.25 
5.34 6.29e-02 44 4 3.00 
6.09 3.18e-02 39 7 3.75 

Km=1; Collision flag=2
4.59 1.11e-01 39 2 2.25 
5.34 6.86e-02 48 4 3.00 
6.09 4.11e-02 36 5 3.75 

Km=2; Collision flag=2
4.59 1.54e-01 54 2 2.25 
5.34 1.89e-01 66 2 3.00 
6.09 9.33e-02 49 3 3.75 

Km=0; Collision flag=1
4.59 5.60e-02 49 5 2.25 
5.34 3.33e-02 35 6 3.00 
6.09 2.06e-02 36 10 3.75 

Km=1; Collision flag=1
4.59 6.86e-02 48 4 2.25 
5.34 4.23e-02 37 5 3.00 
6.09 2.57e-02 36 8 3.75 

Km=2; Collision flag=1
4.59 6.86e-02 36 3 2.25 
5.34 6.29e-02 44 4 3.00 
6.09 4.00e-02 35 5 3.75 

Km=0; Collision flag=0
4.59 3.51e-02 43 7 2.25 
5.34 1.71e-02 36 12 3.00 
6.09 1.54e-02 35 13 3.75 

Km=1; Collision flag=0
4.59 3.33e-02 35 6 2.25 
5.34 2.29e-02 40 10 3.00 
6.09 1.31e-02 39 17 3.75 

Km=2; Collision flag=0
4.59 5.86e-02 41 4 2.25 
5.34 3.71e-02 39 6 3.00 
6.09 2.50e-02 35 8 3.75 

%}

% L(x)=x^2/2+x^3/2 ;R(x)=x^6. Rate=7/12
% BCH outer code (84,98) and LDPC (98,168) base code;
%{

Eb_N0   berr   err sims Es_N0

Km=0; Collision flag=2
6.01 2.63e-01 46 1 3.00 
6.76 7.24e-02 38 3 3.75 
7.51 4.11e-02 36 5 4.50 

Km=1; Collision flag=2
6.01 2.34e-01 41 1 3.00 
6.76 1.03e-01 54 3 3.75 
7.51 5.57e-02 39 4 4.50 

Km=2; Collision flag=2
6.01 2.06e-01 36 1 3.00 
6.76 1.49e-01 52 2 3.75 
7.51 6.29e-02 44 4 4.50 

Km=0; Collision flag=1
6.01 2.14e-01 75 2 3.00 
6.76 4.80e-02 42 5 3.75 
7.51 2.11e-02 37 10 4.50

Km=1; Collision flag=1
6.01 1.09e-01 38 2 3.00 
6.76 3.90e-02 41 6 3.75 
7.51 2.41e-02 38 9 4.50 

Km=2; Collision flag=1
6.01 1.60e-01 56 2 3.00 
6.76 5.71e-02 40 4 3.75 
7.51 3.62e-02 38 6 4.50
 
Km=0; Collision flag=0
6.01 8.76e-02 46 3 3.00 
6.76 1.76e-02 37 12 3.75 

Km=1; Collision flag=0
6.01 2.23e-01 39 1 3.00 
6.76 7.62e-02 40 3 3.75 
7.51 1.67e-02 38 13 4.50 

Km=2; Collision flag=0
6.01 1.86e-01 65 2 3.00 
6.76 7.14e-02 50 4 3.75 
7.51 2.41e-02 38 9 4.50 
%}

% L(x)=x^2/2+x^3/2 ;R(x)=x^8. Rate=11/16=0.6875
% BCH outer code (85,99) and LDPC (99,144) base code;
%{
Km=0; Collision flag=2
Eb_N0   berr   err sims Es_N0
6.49 6.00e-02 42 4 4.25 
6.99 3.02e-02 37 7 4.75 

Km=1; Collision flag=2
6.49 1.17e-01 41 2 4.25 
6.99 5.57e-02 39 4 4.75 

Km=2; Collision flag=2
6.49 1.26e-01 44 2 4.25 
6.99 7.81e-02 41 3 4.75 

Km=0; Collision flag=1
6.49 5.57e-02 39 4 4.25 
6.99 2.41e-02 38 9 4.75 

Km=1; Collision flag=1
6.49 5.14e-02 36 4 4.25 
6.99 3.52e-02 37 6 4.75 

Km=2; Collision flag=1
6.49 9.37e-01 164 1 4.25 
6.99 5.00e-02 35 4 4.75 

Km=0; Collision flag=0
6.49 3.02e-02 37 7 4.25 
6.99 1.00e-02 35 20 4.75 
%}