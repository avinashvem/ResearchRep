#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <stdlib.h>
#include <ctime>

#define num_iter 10  //Number of BP iterations before a forced decimation
#define LLR_MAX 20
#define beta 1.8
#define dc_max 10

#define absi(x) ((x <= 0)? (-1*x):(x))
#define absd(x) ((x <= 0)? (-1.0*x):(x))
#define min(x,y) ((x<=y) ? x:y)
#define max(x,y) ((x>=y) ? x:y)
#define sign(x) ((x <= 0) ? (-1.0):(1.0))

int M, M1, N, N_var, N_chk, tot_edge,L,w,ns;
int dv1,dc1,dv2,dc2;
int dc,dv;

//#include <iostream>
using namespace std;

int parity_check(int c[]);
void read_matrix(FILE* fmat,int no_col,int** fcon);
void message_passing(int** Vcon,int dv,int** Ccon,int dc,float* E_old,float* E,float* LLR,int tot_num_iter);
float diff_sumd(float* vec1,float* vec2,int length);

int main() {
  int *source;
  bool** Mdel;
  int *y_eq;
  float *E,*E_temp,*LLR;
  int *E_hard;
  FILE *Vmat,*Cmat,*param_file;
  int **Vcon, **Ccon;
  float* llr_fin;
  float temp_max;
  int dist_sum,chk_sum;
  float sum;

  int count_n=0,nd,ud;
  int num_err,max_idx;
  int *op;
  int num_chk_err, decim_chki_flag;
  int Mdel_sec_low,Mdel_sec_high;

  int i,j;
  int Dist;

  clock_t init_time=clock();
  clock_t begin_time=clock();
  clock_t end_time;

  param_file=fopen("Simulation_Parameters","r");
  fscanf(param_file,"%d ",&dv1);
  fscanf(param_file,"%d ",&dc1);
  fscanf(param_file,"%d ",&dv2);
  fscanf(param_file,"%d ",&dc2);
  fscanf(param_file,"%d ",&M);
  fscanf(param_file,"%d ",&M1);
  fscanf(param_file,"%d ",&N);
  fscanf(param_file,"%d ",&N_var);
  fscanf(param_file,"%d ",&N_chk);
  fscanf(param_file,"%d ",&tot_edge);
  fscanf(param_file,"%d ",&L);
  fscanf(param_file,"%d ",&w);

  do{
    fscanf(param_file,"\n");
  }while(!feof(param_file));
  
  dc=max(dc1,dc2); dv=dv1+dv2;
  ns=(int) N/(L+w-1);

  printf("M=%d,M1=%d,N=%d,N_var=%d,N_chk=%d,tot_edge=%d,dv=%d,dv1=%d,dv2=%d,dc=%d,dc1=%d,dc2=%d,L=%d,w=%d\n",M,M1,N,N_var,N_chk,tot_edge,dv,dv1,dv2,dc,dc1,dc2,L,w);

  Vcon= new int*[N_var];
  for (i = 0; i < N_var; ++i) {
    Vcon[i] = new int[dv];
  }

  Ccon = new int*[N_chk];
  for (i = 0; i < N_chk; ++i) {
    Ccon[i] = new int[dc];
  }

  Vmat = fopen("Vcon_CC_Coupled_3_6_6_3_200","r");
  read_matrix(Vmat,dv,Vcon);
  Cmat = fopen("Ccon_CC_Coupled_3_6_6_3_200","r");
  read_matrix(Cmat,dc,Ccon);
  
  source=new int[M];
  y_eq=new int[N_var];
  LLR=new float[N_var];
  E=new float[tot_edge+2];
  E_temp=new float[tot_edge+2];
  E_hard=new int[tot_edge+2];
  Mdel=new bool*[N];
  llr_fin=new float[N];
  op=new int[N_var];

  for (i=0;i<N;++i) {
    Mdel[i]=new bool[2];
  }

  count_n=0;
  srand(time(NULL));
  //Generates the source codeword to be Quantized
  for(int i=0;i<M;i++) {
    *(source+i)=(rand()%2);
    *(y_eq+i+N)= *(source+i);
    *(LLR+i+N)=(1-(2*y_eq[i+N]))*beta;
  }
  int Mdel_low=0;
  while(count_n<N) {
    message_passing(Vcon,dv,Ccon,dc,E,E_temp,LLR,num_iter);    

    num_err=0;
    temp_max=0;

    for(i=Mdel_low;i<N;i++){
      if(!Mdel[i][0]){
	Mdel_low=i;
	break;
      }
    }

    Mdel_sec_low=floor(Mdel_low/ns)*ns+1;
    Mdel_sec_high=min(Mdel_sec_low+w*ns,N);
    
    for(i=Mdel_sec_low;i<Mdel_sec_high;++i) {
      if(!Mdel[i][0]) {
    	sum=0;
	for(j=0;j<dv;++j) {
	  sum+=E[Vcon[i][j]];
	}
	*(llr_fin+i)=sum+(*(LLR+i));
	if (temp_max<absd(llr_fin[i])) {
	  max_idx=i;
	  temp_max=absd(llr_fin[i]);
	}
      }
      op[i]=((*(llr_fin+i))<0);
      num_err+=op[i];
    }
    if (temp_max>0) {
      nd=max_idx;
      ud=((rand()/(float)RAND_MAX)>(0.5*(1+tanh(llr_fin[max_idx]))));
    }
    else {
      nd=Mdel_low;
      ud=((rand()/(float)RAND_MAX)>0.5);
    }
    *(LLR+nd)=(1-2*((int)ud))*LLR_MAX;
    Mdel[nd][0]=true;
    Mdel[nd][1]=(bool)ud;
    count_n++;
    for(j=0;j<dv;j++) {
      E_hard[Vcon[nd][j]]=1-2*ud;
    }

    num_chk_err=0;
    // Counting number of checks satisfied
    if ((count_n%100)==0) {
      for(i=0;i<M1;++i) {
	decim_chki_flag=1;	
	chk_sum=0;
	for(j=0;j<dc;++j) {
	  if (Ccon[i][j]!=tot_edge+1) { 
	    if(E_hard[Ccon[i][j]]==0) {
	      decim_chki_flag=0;
	      break;
	    }
	    else {
	      chk_sum+=(int)(0.5*(1-E_hard[Ccon[i][j]]));
	    }
	  }
	  chk_sum=chk_sum%2;
	}
	if(decim_chki_flag) {
	  num_chk_err=num_chk_err+chk_sum;
	}
      }
      printf("Bits decimated/total is %d/%d, ",count_n,N);
      printf("Checks in error/total=%d/%d \n",num_chk_err,M1) ;
    
      Dist=0;
      for(i=M1;i<M1+M;++i) {
	dist_sum=0;
	for(j=0;j<dc1;++j) {
	  if (E_hard[Ccon[i][j]]!=0){
	    dist_sum+=(int)(0.5*(1-E_hard[Ccon[i][j]]));
	  }
	}
	dist_sum=dist_sum%2;
	//		printf("Source bit - %d and Quant bit - %d\n",source[i-M1],sum,abs);
	Dist+=abs(source[i-M1]-dist_sum);
      }
      end_time = clock();
      printf("Distortion is %f \n",(float) Dist/M);
      printf("Time since last 100 decimations is %f\n",float(end_time - begin_time) / CLOCKS_PER_SEC);
      printf("Total Time elapsed is %f\n",float(end_time - init_time) / CLOCKS_PER_SEC);
      begin_time = clock();
    } // End of (count_n%100==0)
  } //End of while loop

  Dist=0;
  for(i=M1;i<M1+M;++i) {
    dist_sum=0;
    for(j=0;j<dc1;++j) {
      dist_sum+=(int)(0.5*(1-E_hard[Ccon[i][j]]));
    }
    dist_sum=dist_sum%2;
    Dist+=abs(source[i-M1]-dist_sum);
  }
  printf("Distortion is %f \n", (float) Dist/M);
}

void message_passing(int** Vcon,int dv,int** Ccon,int dc,float* E_old,float* E,float* LLR,int tot_num_iter) {
  float* temp;
  float f[dc_max];
  float b[dc_max];
  float prod,sum;
  E_old[tot_edge+1]=0;

  int i,j,z;

  for(z=0;z<num_iter;++z) {
    //Bit node Processing
    for (i=0;i<N_var;i++) {
      sum=0;
      for (j=0;j<dv;j++) {
	sum+=E_old[Vcon[i][j]];	
      }
      sum+=LLR[i];
      //printf("sum=%f\n",sum);
      for (j=0;j<dv;j++) {
	if(Vcon[i][j]!=tot_edge+1) {
	  E[Vcon[i][j]]=sum-E_old[Vcon[i][j]];							
	}
	E[Vcon[i][j]]=(E[Vcon[i][j]]>LLR_MAX?LLR_MAX:E[Vcon[i][j]]);
	E[Vcon[i][j]]=(E[Vcon[i][j]]<-LLR_MAX?-LLR_MAX:E[Vcon[i][j]]);
      }
    }
    temp=E_old;
    E_old=E;
    E=temp;

    //Check node Processing
    E_old[tot_edge+1]=LLR_MAX;
   
    for(i=0;i<N_chk;i++) {
      for (j=0;j<dc;j++){
	if(j==0){
	  f[j]=1;
	  b[dc-1-j]=1;
	} else {
	  f[j]=f[j-1]*tanh(E_old[Ccon[i][j-1]]);
	  b[dc-1-j]=b[dc-j]*tanh(E_old[Ccon[i][dc-j]]);
	}
      }
      for (j=0;j<dc;j++){
	E[Ccon[i][j]]=atanh(f[j]*b[j]);
	E[Ccon[i][j]]=(E[Ccon[i][j]]>LLR_MAX?LLR_MAX:E[Ccon[i][j]]);
	E[Ccon[i][j]]=(E[Ccon[i][j]]<-LLR_MAX?-LLR_MAX:E[Ccon[i][j]]);
      }
    }
    temp=E_old;
    E_old=E;
    E=temp;
    E_old[tot_edge+1]=0;
  }
}

float diff_sumd(float* vec1,float *vec2,int length) {
  float sum=0;
  for(int i=0;i<length;++i) {
    sum+=absd(*(vec1+i)-*(vec2+i));
  }
  return(sum);
}

void read_matrix(FILE* fmat,int no_col,int** fcon){
  char tstr[80];
  int l;
  if (fmat == NULL) perror ("Error opening file");
  fgets(tstr,80,fmat);   	
  //printf( tstr );
	
  int node_num=0;
  do{
    for(int i=0;i<no_col;i++) {
      fscanf(fmat,"%d ",&l);	
      *(*(fcon+node_num)+i) = l;
    }
    fscanf(fmat,"\n");
    node_num++;
  } while(!feof(fmat));
  printf("Number of nodes read is %d \n",node_num);
  printf("successfully read the file\n\n");
}
