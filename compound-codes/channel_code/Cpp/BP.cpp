#define N_var 475000//N+M
#define N_chk 237500 //M+M1
#define tot_edge_num 1425000//274500
#define  num_iter 10  //Number of BP iterations before a forced decimation

const int dc= 6;
const int dv= 3;

#define eps 0.078
#define err_resol 1e-3
#define absi(x) ((x <= 0)? (-1*x):(x))
#define absd(x) ((x <= 0)? (-1.0*x):(x))

#define min(x,y) ((x<=y) ? x:y)
#define max(x,y) ((x>=y) ? x:y)
#define sign(x) ((x <= 0) ? (-1.0):(1.0))

#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <stdlib.h>
#include <time.h>
using namespace std;

void read_matrix(FILE* fmat,int no_col,int** fcon);
void message_passing(int** ,int ,int** ,int ,float* ,float* ,int,float );
float diff_sumd(float* vec1,float* vec2,int length);

int main()
{

int *y_eq;
float *LLR;
float* E;
FILE *Vmat,*Cmat;

//int** Vcon; Vcon=(int **)calloc(N_var,sizeof(int *));
int** Vcon = new int*[N_var];
for (int i = 0; i < N_var; ++i)
 {
 	 Vcon[i] = new int[dv];
}

//int** Ccon = new int*[N_chk];
int** Ccon = new int*[N_chk];
for (int i = 0; i < N_chk; ++i)
 {
 	 Ccon[i] = new int[dc];
}

 
Vmat = fopen("Vcon_LDPC.txt","r");
read_matrix(Vmat,dv,Vcon);
Cmat = fopen("Ccon_LDPC.txt","r");
read_matrix(Cmat,dc,Ccon);

float beta=1; // inverse temperature

y_eq=(int*) calloc(N_var,sizeof(int));
LLR=(float*) calloc(N_var,sizeof(float));
E=(float*) calloc(tot_edge_num+1,sizeof(float));

srand(time( NULL ));
for(int i=0;i<N_var;i++) 
	{  
		*(y_eq+i)=((rand()/(float) RAND_MAX)<eps);
		*(LLR+i)=(1-(2*y_eq[i]))*log((1-eps)/eps);
	}
message_passing(Vcon,dv,Ccon,dc,E,LLR,num_iter,beta);	
float* llr_fin;
int op[N_var],ip[N_var],num_err;
llr_fin=(float*)calloc(N_var,sizeof(float));

num_err=0;
int num_err_i=0;
for(int i=0;i<N_var;++i)
 {
	float sum=0;
	for(int j=0;j<dv;++j)
	{
		sum+=E[Vcon[i][j]];
	}
	*(llr_fin+i)=sum+(*(LLR+i));
	op[i]=((*(llr_fin+i))<0);
	ip[i]=((*(LLR+i))<0);
	num_err+=op[i];
	num_err_i+=ip[i];
 }
 printf("Number of errors is %d %d", num_err,num_err_i);
}


void message_passing(int** Vcon,int dv,int** Ccon,int dc,float* E_old,float* LLR,int tot_num_iter,float beta)
{
	float* E; float* temp;
	E=(float*) calloc(tot_edge_num+1,sizeof(float));
	for(int z=0;z<tot_num_iter;++z)
		{
		//Bit node Processing
			
			for (int i=0;i<N_var;++i)
			{
				float sum=0;
				for (int j=0;j<dv;++j)
				{
						sum+=E_old[Vcon[i][j]];	
				}
				sum+=LLR[i];
				for (int k=0;k<dv;++k)
				{
					if(Vcon[i][k]!=tot_edge_num+1)
						{
					float op_ik=sum-E_old[Vcon[i][k]];
					E[Vcon[i][k]]=op_ik;
					//E[Vcon[i][k]]=min(absd(op_ik),25)*sign(op_ik);
						} 
				}
			}
			temp=E_old;
			E_old=E;
			E=temp;
		
		//Check node Processing

			for (int i=0;i<N_chk;++i)
			{
				float prod=1;
				for (int j=0;j<dc;++j)
				{
    				prod*=(tanh(0.5*beta*E_old[Ccon[i][j]]));
				}
				for (int k=0;k<dc;++k)
				{
				if(E_old[Ccon[i][k]]!=0)
					{
						float op_ik=2*atanh(prod/tanh(0.5*E_old[Ccon[i][k]]))/beta;
						E[Ccon[i][k]]=op_ik;
					//	E[Ccon[i][k]]=min(absd(op_ik),25)*sign(op_ik);
					}
					else if(E_old[Ccon[i][k]]==0)							
					{
						E[Ccon[i][k]]=2*atanh(prod)/beta;
					} 
				}
			}
			temp=E_old;
			E_old=E;
   			E=temp;			
//Check the convergence condition: But slowing up quite a bit
			float abs_diff=0;
			if((z%10)==0)
			{ 
				printf("Entered the Convergence check loop \n");
				abs_diff=diff_sumd(E,E_old,tot_edge_num);
				if(abs_diff<err_resol)
					{
						printf("Messages converged...");
						break;
					}
			}
			int op[N_var],num_err=0;
			for(int x=0;x<N_var;++x)
			 {
				 float sum=0;
				for(int j=0;j<dv;++j)
				{
					sum+=E_old[Vcon[x][j]];
				}
				op[x]=(sum+(*(LLR+x)))<0;
				num_err+=op[x];
			 }
			cout << z << ' ' << num_err << endl; 
		}
}


float diff_sumd(float* vec1,float *vec2,int length){
	float sum=0;

		for(int i=0;i<length;++i)
			{
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
			for(int i=0;i<no_col;i++)
			{
				fscanf(fmat,"%d ",&l);	
				*(*(fcon+node_num)+i) = l;
			}
	    	fscanf(fmat,"\n");
		    node_num++;
		}while(!feof(fmat));
	printf("Number of nodes read is %d \n",node_num);	printf("successfully read the file\n\n");
}
