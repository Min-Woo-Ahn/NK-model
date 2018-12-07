#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

const int N=10;
const int K=2;
const int n_k=(int)pow(2,K+1);
const int num_agent=50;
const int Ts=10;
const int Tg=Ts;
const int loop_num=200;
const int Te=1000;

double perf[1000]={0};
double local[1000]={0};

#include "nk_etcfunc.h"
#include "nk_func.h"

main(int argv, char** argc)
{

if(argv!=5)
{
	printf("NK model\n");
	printf("Usage : ./a.out [rho] [lambda] [tau] [option]\n");
	printf("Option : \n");
	printf("1 - Get ensemble average for given condition\n");
	printf("2 - Get time series of organization forms\n");
	printf("3 - Get time series of organization forms with hemming distance\n");
	printf("4 - \n");
	printf("5 - Get performance of individual organization for each time series (local search with Hemming distance)\n");
}

else
{
	int i,j;
	double avr=0,stdev=0,p;
	const int rho=atoi(argc[1]);
	const double lambda=atof(argc[2]);
	const double tau=atof(argc[3]);
	const int option=atoi(argc[4]);
	char filename[200];

	srand(time(NULL));

	if(option==1)
	{
		FILE *f,*g;
		sprintf(filename,"NK-performance-N=%d-K=%d-Ts=%d.csv",N,K,Ts);
		f=fopen(filename,"a");
		sprintf(filename,"rawdata/rawdata-N:%d-K:%d-rho:%d-lamb:%.2lf-tau:%.2lf.csv",N,K,rho,lambda,tau);
		g=fopen(filename,"w");
		fprintf(g,"Simulation condition\n");
		fprintf(g,"\"N=%d, K=%d, (# of agent)=%d, ($ of repetition)=%d\"\n",N,K,num_agent,Te);
		fprintf(g,"\"rho=%d, lambda=%.2lf, tau=%.2lf\"\n\n",rho,lambda,tau);
		fprintf(g,"repetition,firm,average performance during t<=%d\n",loop_num);
		for(i=0;i<Te;i++)
		{
			p=NK_main(rho,lambda,tau);
			avr+=p;
			stdev+=p*p;
			for(j=0;j<num_agent;j++)
			{
				fprintf(g,"%d,%d,%lf\n",i+1,j,perf[j]);
			}
		}
		avr/=Te;
		stdev=sqrt(stdev/Te-avr*avr);
		printf("Average fitness for rho=%d, lambda=%lf, tau=%lf : %lf\n",rho,lambda,tau,avr);
		fprintf(f,"%d,%lf,%lf,%lf,%lf,%lf,%lf\n",rho,lambda,tau,avr,stdev,avr-2*stdev/Te,avr+2*stdev/Te);
		fclose(f);
		fclose(g);
	}

	else if(option==2)
	{
		FILE *f,*g;
		sprintf(filename,"NK-hemming-performance-N=%d-K=%d.csv",N,K);
		f=fopen(filename,"a");
		sprintf(filename,"rawdata/rawdata-hemming-N:%d-K:%d-rho:%d-lamb:%.2lf-tau:%.2lf.csv",N,K,rho,lambda,tau);
		g=fopen(filename,"w");
		fprintf(g,"Simulation condition\n");
		fprintf(g,"\"N=%d, K=%d, (# of agent)=%d, ($ of repetition)=%d\"\n",N,K,num_agent,Te);
		fprintf(g,"\"rho=%d, lambda=%.2lf, tau=%.2lf\"\n\n",rho,lambda,tau);
		fprintf(g,"repetition,firm,average performance during t<=%d\n",loop_num);
		for(i=0;i<Te;i++)
		{
			p=NK_hemming(rho,lambda,tau,1);
			avr+=p;
			stdev+=p*p;
			for(j=0;j<num_agent;j++)
			{
				fprintf(g,"%d,%d,%lf\n",i+1,j,perf[j]);
			}
		}
		avr/=Te;
		stdev=sqrt(stdev/Te-avr*avr);
		printf("Average fitness for rho=%d, lambda=%lf, tau=%lf : %lf\n",rho,lambda,tau,avr);
		fprintf(f,"%d,%lf,%lf,%lf,%lf,%lf,%lf\n",rho,lambda,tau,avr,stdev,avr-2*stdev/Te,avr+2*stdev/Te);
		fclose(f);
		fclose(g);
	}
	else if(option==3)
	{
		FILE *g;
		sprintf(filename,"data_emergent/rawdata-hemming-time-N:%d-K:%d-rho:%d-lamb:%.2lf-tau:%.2lf.csv",N,K,rho,lambda,tau);
		g=fopen(filename,"w");
		fprintf(g,"Simulation condition\n");
		fprintf(g,"\"N=%d, K=%d, (# of agent)=%d, ($ of repetition)=%d\"\n",N,K,num_agent,Te);
		fprintf(g,"\"rho=%d, lambda=%.2lf, tau=%.2lf\"\n\n",rho,lambda,tau);
		fprintf(g,"timestep,performance\n");
		for(i=0;i<Te;i++)
		{
			for(j=0;j<loop_num;j++) perf[j]=0;
			NK_hemming(rho,lambda,tau,3);
			for(j=0;j<loop_num;j++)
			{
				fprintf(g,"%d,%lf\n",j,perf[j]);
			}
		}
		fclose(g);

	}

	else if(option==20)
	{
		int T=500;
		double l[3]={0};
		FILE *f;
		for(i=0;i<T;i++)
		{
			NK_hemming_local(rho,lambda,tau);
			for(j=0;j<3;j++)
			{
				l[j]+=local[j]/(double)T;
			}
		}
		sprintf(filename,"data/localwalk-N:%d-K:%d.csv",N,K);
		f=fopen(filename,"a");
		fprintf(f,"%d,%lf,%lf,%lf,%lf,%lf\n",rho,lambda,tau,local[2],local[1],local[0]);
		fclose(f);
	}
	else if(option==4)
	{
		FILE *f;
		double n_local=0,n_hemming=0,f_local=0,f_hemming=0;
		for(i=0;i<500;i++)
		{
			NK_landscape_hemming(rho,lambda,tau);
			n_local+=perf[0]/500;
			n_hemming+=perf[1]/500;
			f_local+=perf[2]/500;
			f_hemming+=perf[3]/500;
		}
		sprintf(filename,"localopt-fitness-N:%d-K:%d.csv",N,K);
		f=fopen(filename,"a");
		fprintf(f,"%d,%lf,%lf,%lf,%lf,%lf,%lf\n",rho,lambda,tau,n_local,n_hemming,f_local,f_hemming);
		fclose(f);
	}
	else if(option==5) NK_hemming(rho,lambda,tau,2);
	else if(option==6) NK_sd_fitness(rho,lambda,tau,2);
	else if(option==7) NK_sd_fitness(rho,lambda,tau,3);
	else if(option==8)
	{
		double same=0;
		for(i=0;i<Te;i++) same+=NK_sd_fitness(rho,lambda,tau,4);
		printf("average # number of global optimum = %lf\n",same/Te);
	}
	else if(option==9)
	{
		FILE *f;
		double p;
		double p_global=0;
		double p_stdev=0;
		for(i=0;i<Te;i++)
		{
			p=NK_count_global_optimum(rho,lambda,tau);
			p_global+=p/Te;
			p_stdev+=p*p/Te;
		}
		p_stdev=sqrt(p_stdev-p_global*p_global);
		sprintf(filename,"global_optimum_p-N:%d-K:%d.csv",N,K);
		f=fopen(filename,"a");
		fprintf(f,"%d,%lf,%lf,%lf,%lf\n",rho,lambda,tau,p_global,p_stdev);
		fclose(f);
	}
}
}
