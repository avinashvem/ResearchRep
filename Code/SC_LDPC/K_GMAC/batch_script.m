clear all
addpath('../PEG/')
addpath('~/Github/avinash/ResearchRep/Code/SC_LDPC/for_Avinash/');
chcode.N=28500;
chcode.M=144;
chcode.kmsg=84;
chcode.tBCH=0;
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

params.EsN0_vec=[2.75];
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
                errMat=main_extra_missed_perms(Kparams,params,chcode);
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
%-----------------Results-------------------
%Kb=10
%{
Ka=100; Kb=10 (110 total);Km=0; Collision flag=0
Eb_N0   berr   err sims Es_N0
1.01    1.0    100  1    -2
1.51    0.91   91  1  
2.01    0.43   43  1  
2.51    0.07   14  2  
3.01  1.83e-2  11  6
3.51  7.00e-3   7 10
4.01  3.00e-3   3 10
4.51  4.00e-3   4 10
5.01  3.00e-3   3 10

Ka=100; Kb=10 (110 total);Km=1; Collision flag=0
Eb_N0   berr   err sims Es_N0
1.01    0.99    100  1    -2
1.51    0.86   91  1  
2.01    0.25   43  1  
2.51  9.50e-2  19  2  
3.01  2.80e-2  14  5
3.51  2.50e-2  10  4
4.01  1.00e-2   4  4

Ka=100; Kb=10 (110 total);Km=2; Collision flag=0
Eb_N0   berr err sims Es_N0
1.01 9.90e-01 99 1 -2.00 
1.51 9.00e-01 90 1 -1.50 
2.01 5.60e-01 56 1 -1.00 
2.51 7.00e-02 14 2 -0.50 
3.01 5.50e-02 11 2 0.00 
3.51 2.75e-02 11 4 0.50 
4.01 2.40e-02 12 5 1.00 
4.51 2.20e-02 11 5 1.50 


Ka=100; Kb=10 (110 total);Km=0; Collision flag=1
Eb_N0   berr   err sims Es_N0
1.01    0.99   99  1    -2
1.51    0.89   89  1  
2.01    0.47   47  1  
2.51    0.12   12  1  
3.01  5.00e-2  10  2
3.51  2.75e-2  11  4
4.01  3.00e-2  12  4
4.51  2.20e-2  11  5

Ka=100; Kb=10 (110 total);Km=1; Collision flag=1
colliding index (index is Ka/2=50) is input wrong to BP decoder 
Eb_N0   berr  err sims Es_N0
1.01 1.00e+00 100 1 -2.00 
1.51 9.10e-01 91 1 -1.50 
2.01 3.00e-01 30 1 -1.00 
2.51 1.50e-01 15 1 -0.50 
3.01 5.00e-02 10 2 0.00 
3.51 3.00e-02 15 5 0.50 
4.01 3.25e-02 13 4 1.00 
4.51 2.50e-02 10 4 1.50 

Ka=100; Kb=10 (110 total);Km=1; Collision flag=1
colliding indices (twice) is input correct to BP decoder. Random index is input wrong
Eb_N0  berr  err sims Es_N0
1.01 1.00e+00 100 1 -2.00 
1.51 9.40e-01 94 1 -1.50 
2.01 6.70e-01 67 1 -1.00 
2.51 1.90e-01 19 1 -0.50 
3.01 6.00e-02 12 2  0.00 
3.51 4.00e-02 12 3  0.50 
4.01 3.33e-02 10 3  1.00 

Ka=100; Kb=10 (110 total);Km=2; Collision flag=1
Random index is input wrong. colliding index (twice) is input correct to BP decoder.
Eb_N0  berr  er sims Es_N0
1.01 9.80e-01 98 1 -2.00 
1.51 9.60e-01 96 1 -1.50 
2.01 3.70e-01 37 1 -1.00 
2.51 1.80e-01 18 1 -0.50 
3.01 7.50e-02 15 2  0.00 
3.51 5.00e-02 15 3  0.50 
4.01 4.33e-02 13 3  1.00
%}

% Ka=175Kb=0.1*Ka; Various rates
%{
L(x)=x^2/2+x^3/2 ;R(x)=x^6
Rate=7/12=0.5833
(84,144) base code; Km=0; Collision flag=0
Eb_N0   berr err sims Es_N0
4.34 6.14e-02 43 4 2.00 
5.09 2.86e-02 35 7 2.75 
5.84 1.14e-02 36 18 3.50 
6.59 7.43e-03 26 20 4.25 

Ka=175,Kb=10
(84,144) base code; Km=1; Collision flag=0
4.34 5.29e-02 37 4 2.00 
5.09 3.27e-02 40 7 2.75 
5.84 2.41e-02 38 9 3.50 
6.59 1.32e-02 37 16 4.25 

Ka=175,Kb=10
(84,144) base code; Km=0; Collision flag=2
4.34 7.81e-02 41 3 2.00 
5.09 3.52e-02 37 6 2.75 
5.84 3.43e-02 36 6 3.50 
6.59 1.58e-02 36 13 4.25
 
Ran with maxIters=1000 (to solve num_err=1,2,.., issue). Didn't seem to help
4.59 1.26e-01 44 2 2.25 
5.34 6.67e-02 35 3 3.00 
6.09 3.71e-02 39 6 3.75 

(84,144) base code; Km=1; Collision flag=2
4.34 1.69e-01 59 2 2.00 
5.09 8.19e-02 43 3 2.75 
5.84 4.00e-02 42 6 3.50 
6.59 2.71e-02 38 8 4.25 

(84,144) base code; Km=2; Collision flag=2
4.34 1.06e-01 37 2 2.00 
5.09 6.86e-02 36 3 2.75 
5.84 6.86e-02 36 3 3.50 
6.59 3.10e-02 38 7 4.25 

(84,144) base code; Km=0; Collision flag=1
4.34 9.71e-02 51 3 2.00 
5.09 6.29e-02 44 4 2.75 
5.84 3.18e-02 39 7 3.50 
6.59 2.57e-02 36 8 4.25 

(84,144) base code; Km=1; Collision flag=1
4.34 6.86e-02 36 3 2.00 
5.09 5.00e-02 35 4 2.75 
5.84 3.18e-02 39 7 3.50 
6.59 2.57e-02 36 8 4.25 


L(x)=x^2/2+x^3/2 ;R(x)=x^7
(90,140) base code; Km=0; Collision flag=0
Rate=9/14=0.6429
Eb_N0   berr   err sims Es_N0
4.92  5.42e-2   38  4    3.00
5.67  1.87e-2   36 11    3.75
6.42  1.11e-2   37 19    4.50


L(x)=x^2/2+x^3/2 ;R(x)=x^8
(88,128) base code; Km=0; Collision flag=0
Rate=11/16=0.6875
Eb_N0   berr   err sims Es_N0
3.63    1.0    172  1    2.00
4.38    0.13    70  3    2.75
5.13  7.05e-2   37  3    3.50
5.88  3.02e-2   37  7    4.25  


L(x)=x^2/2+x^3/2 ;R(x)=x^9
(91,126) base code; Km=0; Collision flag=0
Rate=13/18=0.7222
Eb_N0   berr   err sims Es_N0
3.41    1.0    174  1    2.00
4.16    0.98   171  1    2.75
4.91    0.98   171  1    3.50
5.66  5.57e-2   39  4    4.25  

%}
%--------------Various Kb--------------------
%{
Ka=100; Kb=3 (103 total); Collision flag=0
Eb_N0   berr   err sims Es_N0
1.01    0.98   98  1    -2
1.51    0.84   84  1  
2.01    0.27   27  1  
2.51    0.07   14  2  
3.01  1.71e-2  12  7
3.51  2.00e-3   2 10
4.01    0       0 10

Ka=100; Kb=5 (105 total); Collision flag=0
Eb_N0   berr   err sims Es_N0
1.01    0.99   99  1    -2
1.51    0.78   78  1  
2.01    0.18   18  1  
2.51    0.06   12  2  
3.01  2.20e-2  11  5
3.51  4.00e-3   4 10
4.01    0       0 10


Ka=100; Kb=20 (120 total); Collision flag=0
Eb_N0   berr   err sims Es_N0
1.01    1.0    100  1    -2
1.51    0.96   91   1  
2.01    0.36   43   1  
2.51    0.14   14   2  
3.01  4.00e-2  12   3
3.51  1.57e-3  11   7
4.01  1.25e-2  10   8
4.51  7.00e-3   7  10
5.01  2.00e-3   2 10


%---------With one collision----------------------------------------------------
Ka=100; Kb=3 (103 total); Collision flag=1
Eb_N0   berr   err sims Es_N0
1.01    0.96   96  1    -2
1.51    0.82   82  1  
2.01    0.26   26  1  
2.51    0.08   16  2  
3.01  2.75e-2  11  4
3.51  2.00e-2  10  5
4.01  2.00e-2  10  5
4.51  2.00e-2  10  5

Ka=100; Kb=5 (105 total); Collision flag=1
Eb_N0   berr   err sims Es_N0
1.01    0.97   97  1    -2
1.51    0.85   85  1  
2.01    0.26   26  1  
2.51    0.10   10  1  
3.01  3.66e-2  11  3
3.51  2.00e-2  10  5
4.01  2.00e-2  10  5
4.51  2.00e-2  10  5

Ka=100; Kb=20 (120 total); Collision flag=0
Eb_N0   berr   err sims Es_N0
1.01    1.00  100  1    -2
1.51    0.97   97  1  
2.01    0.74   74  1  
2.51    0.24   24  1  
3.01  5.50e-2  11  2
3.51  4.00e-2  12  3
4.01  2.75e-2  11  4
4.51  4.66e-2  14  3
%}
%{
results.Ka=Kparams.Ka;
    results.Kb=Kparams.Kb;
    results.Km=Kparams.Km;
    results.collision=params.collision_flag;
    results.EsN0_vec=params.EsN0_vec;
    results.N=chcode.N;
    results.M=chcode.M;
    results.helpString='The rows represetn Es/N0 values. Eb/N0 =s Es/N0+3.01dB. The 8 columns correspond respectively to EbN0 bitErrPrb number_errors num_simulations Es/N0 Ka Kb Km';
    results.data=errMat;
    load('chcode_results.mat')
    rLength=length(resultsMat);
    resultsMat(rLength+1)=results;
    save('chcode_results.mat','resultsMat');
%}