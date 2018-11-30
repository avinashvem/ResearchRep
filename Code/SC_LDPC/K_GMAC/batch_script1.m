clear all
addpath('../PEG/')
addpath('~/Github/avinash/ResearchRep/Code/SC_LDPC/for_Avinash/');
chcode.N=28500;
chcode.M=172;
chcode.L=[0.5 0.5;2 3];chcode.R=[1;5];
chcode.base_rate=1/2;
chcode.repeat_flag=0;
chcode.repeat=1;
girth=0;
%--Michael's PEG construction------------------
while girth<8
[chcode.Vcon,chcode.Ccon,chcode.eMax,girth]=LDPC_PEG_cnstr(chcode.M,round(chcode.M*(1-chcode.base_rate)),chcode.L,chcode.R,chcode.repeat_flag);
end

Ka_vec=[50 75];
Km_vec=[0 1 2];
collision_vec=[0 1 2];

params.EsN0_vec=sort(-2:0.5:1);
params.missing_idx=[];

op=zeros(length(Ka_vec),length(Km_vec),length(collision_vec),length(params.EsN0_vec),8);
for idxK=1:length(Ka_vec)
    Kparams.Ka=Ka_vec(idxK);
    Kparams.Kb=ceil(Kparams.Ka/10);
    
    for idxc=1:length(collision_vec)
        params.collision_flag=collision_vec(idxc);
        
        for idxm=1:length(Km_vec)
            Kparams.Km=Km_vec(idxm);
            load('chcode_results.mat')
            keyvalue=[num2str(Kparams.Ka) ',' num2str(Kparams.Kb) ',' num2str(Kparams.Km) ',' num2str(params.collision_flag)];
 
            prverrMat=[];
            
            if Hash.isKey(keyvalue)
                fprintf('Already partially computed.\n');
                prverrMat=Hash(keyvalue);
                
                if max(prverrMat(:,5))>=min(params.EsN0_vec)
                    fprintf('Check the input SNR values. stored errMat for key %s is',keyvalue);
                    Hash(keyvalue)
                    continue
                end
            end

            fprintf('Parameters: Ka=%d,Kb=%d,Km=%d,collision flag=%d\n',Kparams.Ka,Kparams.Kb,Kparams.Km,params.collision_flag)
            errMat=main_extra_missed_perms(Kparams,params,chcode);

            if numel(prverrMat)>0
                errMat=[prverrMat;errMat]
            end

            load('chcode_results.mat')
            Hash(keyvalue)=errMat;
            assert(aN==chcode.N);
            assert(aM==chcode.M);
            save('chcode_results.mat','Hash','aN','aM')

            for i=1:numel(errMat(:,1))
                fprintf('%3.2f %3.2e %d %d %3.2f \n',errMat(i,1),errMat(i,2),errMat(i,3),errMat(i,4),errMat(i,5));
            end
            
        end
      end
end