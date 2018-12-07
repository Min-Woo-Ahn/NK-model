double calculate_cont(int *seq, double** c, int** inter)
{
	double cont=0;
	for(int i=0;i<N;i++) cont+=c[i][get_number_from_binary(i,seq,inter)];
	cont/=N;
	return cont;
}

void calculate_fitness(double* fitness,double** c, int** inter)
{
	int *binary;
	int n=(int)pow(2.0,N);
	binary=new int[N];
	for(int i=0;i<n;i++)
	{
		get_binary_from_number(i,binary);
		fitness[i]=calculate_cont(binary,c,inter);
	}
}

int distance(int* x,int* y)
	{
	int dist=0,i;
	for(i=0;i<N;i++)
	{
		if(x[i]!=y[i])dist++;
	}
	return dist;
}

int local_optimum(int *seq, double **c, int** inter)
{
        int i,*temp;
        double fitness;
        temp=new int[N];
        for(i=0;i<N;i++)temp[i]=seq[i];
        fitness=calculate_cont(temp,c,inter);
        for(i=0;i<N;i++)
	{
                temp[i]=1-temp[i];
                if(fitness<calculate_cont(temp,c,inter))return 0;
                temp[i]=1-temp[i];
        }
        return 1;
}

int local_optimum(int *seq, double **c, int** inter,int *global)
{
        int i,*temp;
        double fitness;
        temp=new int[N];
        for(i=0;i<N;i++) temp[i]=seq[i];
        fitness = calculate_cont(temp,c,inter) * hemming_distance(seq,global,N);
        for(i=0;i<N;i++)
	{
                temp[i]=1-temp[i];
                if(fitness<calculate_cont(temp,c,inter)*hemming_distance(temp,global,N))return 0;
                temp[i]=1-temp[i];
        }
        return 1;
}

class agent
{
	public:
		int *seq;
		double cont;
	agent()
	{
		cont=0;
	}
	void sequence()
	{
		int i;
		seq=new int[N];
		for(i=0;i<N;i++) seq[i]=rand()%2;
	}
	void calculate_contribution(double** c,int** inter)
	{
		cont=calculate_cont(seq,c,inter);
	}
	void local_search(double** c, int** inter, int* fix, int rho)
	{
		int x;
		do
		{
			x=rand()%N;
		}while(fixed_decision(x,fix,rho));
		calculate_contribution(c,inter);
		double prev=cont;
		seq[x]=1-seq[x];
		calculate_contribution(c,inter);
		if(prev>=cont)
		{
			seq[x]=1-seq[x];
			cont=prev;
		}
	}
	void local_search2(double** c, int** inter)
	{
		int i;
		int highest_i=-1;
		calculate_contribution(c,inter);
		double highest_cont=cont;
		for(i=0;i<N;i++)
		{
			seq[i]=1-seq[i];
			calculate_contribution(c,inter);
			if(highest_cont<cont)
			{
				highest_cont=cont;
				highest_i=i;
				seq[i]=1-seq[i];
			}
			else
			{
				seq[i]=1-seq[i];
			}
		}
		if(highest_i!=-1)seq[highest_i]=1-seq[highest_i];
		calculate_contribution(c,inter);
	}
	int local_search_hemming(double** c, int* global_future, int** inter,int *fix, int rho)
	{
		int x;
		double hf_prev,hf;
		do	
		{
			x=rand()%N;
		}while(fixed_decision(x,fix,rho));
		calculate_contribution(c,inter);
		double prev=cont;
		hf_prev=cont*hemming_distance(seq,global_future,N)/N;

		seq[x]=1-seq[x];
		calculate_contribution(c,inter);
		hf=cont*hemming_distance(seq,global_future,N)/N;
		if(hf_prev>=hf)
		{
			seq[x]=1-seq[x];
			cont=prev;
			return -1;
		}
		if(cont<prev) return 1;
		else return 0;
	}
	void local_random(double** c,int** inter)
	{
		int i=rand()%N;
		seq[i]=1-seq[i];
		calculate_contribution(c,inter);
	}
	~agent()
	{
		delete[] seq;
        }
};

void get_global_optimum(int* seq, double** c, int** inter)
{
	double cont,cont_max=0;
	int i,num_max=0;
	const int num=(int)pow(2.0,N);
	int *config;
	config=new int[N];
	for(i=0;i<num;i++)
	{
		get_binary_from_number(i,config);
		cont=calculate_cont(config,c,inter);
		if(cont>cont_max)
		{
			cont_max=cont;
			num_max=i;
		}
	}
	get_binary_from_number(num_max,seq);
	delete[] config;
}

void interaction_allocation(int **inter)
{
        int i,j;
        for(i=0;i<N;i++)
	{
                for(j=0;j<K;j++)
		{
                        do
			{
                                inter[i][j]=rand()%N;
                        }while(redundancy(inter[i],i,j+1));
                }
        }
}

void reallocate_fitness(double **c,double tau)
{
        int i,j;
        for(i=0;i<N;i++)
	{
                for(j=0;j<n_k;j++)
		{
			c[i][j] = c[i][j] * (1-tau) + tau * uniform();
                }
        }
}

double NK_count_global_optimum(int rho,double lambda,double tau)
{
	int i,j;
	double **c;
	int **interaction;
	int global_count=0;
	agent *a;
	a=new agent[num_agent];
	for(i=0;i<num_agent;i++)a[i].sequence();

	c=new double*[N];
	for(i=0;i<N;i++)
	{
		c[i]=new double[n_k];
		for(j=0;j<n_k;j++)c[i][j]=uniform();
	}

	interaction=new int*[N];
	for(i=0;i<N;i++)interaction[i]=new int[K];
	interaction_allocation(interaction);

	int *global;
	global=new int[N];

	int* list;
	list=new int[N];
	int *temp;	
	temp=new int[N];
	for(i=0;i<N;i++)list[i]=i;
	for(i=0;i<num_agent;i++) perf[i]=0;
	for(int loop=0;loop<loop_num;loop++)
	{
		if(loop%Ts==0&&tau!=0)
		{
			reallocate_fitness(c,tau);
		}
		if((loop%Tg==0&&tau!=0)||loop==0)
		{
			get_global_optimum(global,c,interaction);
			shuffle(list);
			for(i=0;i<num_agent;i++)
			{
				for(j=0;j<rho;j++)
				{
					if(uniform()<=lambda) a[i].seq[list[j]]=global[list[j]];
					else a[i].seq[list[j]]=1-global[list[j]];
				}
			}
		}
		for(i=0;i<num_agent;i++) a[i].local_search(c,interaction,list,rho);
		for(i=0;i<num_agent;i++) a[i].calculate_contribution(c,interaction);

		if(loop%Ts==(Ts-1))
		{
			for(i=0;i<num_agent;i++)
			{
				if(distance(global,a[i].seq)==0)global_count++;
			}
		}
	}
	for(i=0;i<N;i++)
	{
		delete[] c[i];
		delete[] interaction[i];
	}
	delete[] c;
	delete[] interaction;
	delete[] a;
	delete[] global;
	delete[] list;
	return (double)global_count/(num_agent*(loop_num/Ts));
}

double NK_main(int rho,double lambda,double tau)
{
	int i,j;
	double **c;
	int **interaction;
	agent *a;
	a=new agent[num_agent];
	for(i=0;i<num_agent;i++) a[i].sequence();

	c=new double*[N];
	for(i=0;i<N;i++)
	{
		c[i]=new double[n_k];
		for(j=0;j<n_k;j++) c[i][j]=uniform();
	}

	interaction=new int*[N];
	for(i=0;i<N;i++)interaction[i]=new int[K];
	interaction_allocation(interaction);

	int *global;
	global=new int[N];

	int* list;
	list=new int[N];
	int *temp;	
	temp=new int[N];
	for(i=0;i<N;i++)list[i]=i;
	for(i=0;i<num_agent;i++) perf[i]=0;
	for(int loop=0;loop<loop_num;loop++)
	{
		if(loop%Ts==0&&tau!=0)
		{
			reallocate_fitness(c,tau);
		}
		if((loop%Tg==0&&tau!=0)||loop==0)
		{
			get_global_optimum(global,c,interaction);
			shuffle(list);
			for(i=0;i<num_agent;i++)
			{
				for(j=0;j<rho;j++)
				{
					if(uniform()<=lambda) a[i].seq[list[j]]=global[list[j]];
					else a[i].seq[list[j]]=1-global[list[j]];
				}
			}
		}
		for(i=0;i<num_agent;i++) a[i].local_search(c,interaction,list,rho);
		for(i=0;i<num_agent;i++) a[i].calculate_contribution(c,interaction);
		for(i=0;i<num_agent;i++) perf[i]+=a[i].cont/loop_num;
	}
	double avr=0;
	for(i=0;i<num_agent;i++)
	{
		avr+=a[i].cont;
	}
	for(i=0;i<N;i++)
	{
		delete[] c[i];
		delete[] interaction[i];
	}
	delete[] c;
	delete[] interaction;
	delete[] a;
	delete[] global;
	delete[] list;
	return avr/num_agent;
}


double NK_hemming(int rho, double lambda, double tau, int option)
{
	int i,j;
	double **c,**c_future;
	int **interaction;
	FILE *f;
	char filename[200];
	agent *a;
	a=new agent[num_agent];
	for(i=0;i<num_agent;i++) a[i].sequence();

	for(i=0;i<num_agent;i++) perf[i]=0;

	c=new double*[N];
	c_future=new double*[N];
	for(i=0;i<N;i++)
	{
		c[i]=new double[n_k];
		c_future[i]=new double[n_k];
		for(j=0;j<n_k;j++) c[i][j]=uniform();
		for(j=0;j<n_k;j++) c_future[i][j]=c[i][j]*(1-tau)+tau*uniform();
	}

	interaction=new int*[N];
	for(i=0;i<N;i++) interaction[i]=new int[K];
	interaction_allocation(interaction);

	int *global,*global_future;
	global=new int[N];
	global_future=new int[N];
	get_global_optimum(global_future,c_future,interaction);

	int* list;
	list=new int[N];
	int *temp;	
	temp=new int[N];
	for(i=0;i<N;i++)list[i]=i;
	if(option==2)
	{
		sprintf(filename,"data/hemming-series-N:%d-K:%d-%d-%.2lf-%.2lf.csv",N,K,rho,lambda,tau);
		f=fopen(filename,"w");
		fprintf(f,"Firm ID");
		for(i=0;i<num_agent;i++) fprintf(f,",%d",i);
		fprintf(f,"\n");
	}

	for(int loop=0;loop<loop_num;loop++)
	{
		if(loop%Ts==0&&tau!=0)
		{
			for(i=0;i<N;i++)
			{
				for(j=0;j<n_k;j++)
				{
					c[i][j]=c_future[i][j];
					c_future[i][j]=c[i][j]*(1-tau)+tau*uniform();
				}
			}
			get_global_optimum(global_future,c_future,interaction);
		}
		if((loop%Tg==0&&tau!=0)||loop==0)
		{
			get_global_optimum(global,c,interaction);
			shuffle(list);
			for(i=0;i<num_agent;i++)
			{
				for(j=0;j<rho;j++)
				{
					if(uniform()<=lambda) a[i].seq[list[j]]=global[list[j]];
					else a[i].seq[list[j]]=1-global[list[j]];
				}
			}
		}
		for(i=0;i<num_agent;i++) a[i].local_search_hemming(c,global_future,interaction,list,rho);
		for(i=0;i<num_agent;i++) a[i].calculate_contribution(c,interaction);

		if(option==1)
		{
			for(i=0;i<num_agent;i++) perf[i]+=a[i].cont/loop_num;
		}
		else if(option==2)
		{
			fprintf(f,"T=%d",loop);
			for(i=0;i<num_agent;i++) fprintf(f,",%d",(int)hemming_distance(a[i].seq,global_future,N));
			fprintf(f,"\n");
		}
		else if(option==3)
		{
			for(i=0;i<num_agent;i++) perf[loop]+=a[i].cont/num_agent;
		}
	}
	if(option==2) fclose(f);
	double avr=0;
	for(i=0;i<num_agent;i++) avr+=perf[i];

	for(i=0;i<N;i++)
	{
		delete[] c[i];
		delete[] interaction[i];
	}
	delete[] c;
	delete[] interaction;
	delete[] a;
	delete[] global;
	delete[] global_future;
	delete[] list;
	return avr/num_agent;
}

double NK_hemming_local(int rho,double lambda,double tau)
{
	int i,j;
	double **c,**c_future;
	int **interaction;
	FILE *f;
	char filename[200];
	agent *a;
	a=new agent[num_agent];
	for(i=0;i<num_agent;i++) a[i].sequence();

	for(i=0;i<num_agent;i++) perf[i]=0;

	c=new double*[N];
	c_future=new double*[N];
	for(i=0;i<N;i++)
	{
		c[i]=new double[n_k];
		c_future[i]=new double[n_k];
		for(j=0;j<n_k;j++)c[i][j]=uniform();
		for(j=0;j<n_k;j++)c_future[i][j]=c[i][j]*(1-tau)+tau*uniform();
	}

	interaction=new int*[N];
	for(i=0;i<N;i++) interaction[i]=new int[K];
	interaction_allocation(interaction);

	int *global,*global_future;
	global=new int[N];
	global_future=new int[N];
	get_global_optimum(global_future,c_future,interaction);

	int* list;
	list=new int[N];
	int *temp;	
	temp=new int[N];
	for(i=0;i<N;i++)list[i]=i;
	for(int loop=0;loop<loop_num;loop++)
	{
		if(loop%Ts==0&&tau!=0)	
		{
			for(i=0;i<N;i++)
			{
				for(j=0;j<n_k;j++)
				{
					c[i][j]=c_future[i][j];
					c_future[i][j]=c[i][j]*(1-tau)+tau*uniform();
				}
			}
			get_global_optimum(global_future,c_future,interaction);
		}
		if((loop%Tg==0&&tau!=0)||loop==0)
		{
			get_global_optimum(global,c,interaction);
			shuffle(list);
			for(i=0;i<num_agent;i++)
			{
				for(j=0;j<rho;j++)
				{
					if(uniform()<=lambda) a[i].seq[list[j]]=global[list[j]];
					else a[i].seq[list[j]]=1-global[list[j]];
				}
			}
		}
		for(i=0;i<num_agent;i++) a[i].local_search_hemming(c,global_future,interaction,list,rho);
		for(i=0;i<num_agent;i++) a[i].calculate_contribution(c,interaction);
		for(i=0;i<num_agent;i++) local[i]+=hemming_distance(a[i].seq,global_future,N)/(loop_num*N);
	}
	double avr=0;
	for(i=0;i<num_agent;i++) avr+=perf[i];

	for(i=0;i<N;i++)
	{
		delete[] c[i];
		delete[] interaction[i];
	}
	delete[] c;
	delete[] interaction;
	delete[] a;
	delete[] global;
	delete[] global_future;
	delete[] list;
	return avr/num_agent;
}

void NK_landscape_hemming(int rho,double lambda,double tau)
{
	int i,j;
	double **c,**c_future;
	int **interaction;
	char fn1[200],fn2[200];

	c=new double*[N];
	c_future=new double*[N];
	for(i=0;i<N;i++)
	{
		c[i]=new double[n_k];
		c_future[i]=new double[n_k];
		for(j=0;j<n_k;j++)c[i][j]=uniform();
		for(j=0;j<n_k;j++)c_future[i][j]=c[i][j]*(1-tau)+tau*uniform();
	}
	for(i=0;i<4;i++) perf[i]=0;
	interaction=new int*[N];
	for(i=0;i<N;i++) interaction[i]=new int[K];
	interaction_allocation(interaction);

	int *global,*global_future;
	global=new int[N];
	global_future=new int[N];

	int *temp;	
	temp=new int[N];
	int pow2=(int)pow(2,N);
	double *local_opt,*local_hemming;
	local_opt=new double[pow2/10];
	local_hemming=new double[pow2/10];
	int *d_local,*d_hemming;
	d_local=new int[pow2/10];
	d_hemming=new int[pow2/10];
	int num_local,num_hemming;
	FILE *f,*g;
	sprintf(fn1,"data/local-optimum-N:%d-K:%d-%d-%.2lf-%.2lf.csv",N,K,rho,lambda,tau);
	sprintf(fn2,"data/local-optimum-hamming-N:%d-K:%d-%d-%.2lf-%.2lf.csv",N,K,rho,lambda,tau);
	f=fopen(fn1,"w");
	g=fopen(fn2,"w");
	for(int loop=0;loop<loop_num;loop++)
	{
		if(loop%Ts==0)
		{
			for(i=0;i<N;i++)
			{
				for(j=0;j<n_k;j++)
				{
					c[i][j]=c_future[i][j];
					c_future[i][j]=c[i][j]*(1-tau)+tau*uniform();
				}
			}
			get_global_optimum(global_future,c_future,interaction);
			get_global_optimum(global,c,interaction);
			num_local=0;
			num_hemming=0;
			for(i=0;i<pow2;i++)
			{
				get_binary_from_number(i,temp);
				if(local_optimum(temp,c,interaction,global))
				{
					local_hemming[num_hemming]=calculate_cont(temp,c,interaction);
					d_hemming[num_hemming]=distance(temp,global);
					num_hemming++;
				}
				if(local_optimum(temp,c,interaction))
				{
					local_opt[num_local]=calculate_cont(temp,c,interaction);
					d_local[num_local]=distance(temp,global);
					num_local++;
				}
			}
			double ratio=loop_num/Ts;
			perf[0]+=(double)num_local/ratio;
			perf[1]+=(double)num_hemming/ratio;
			perf[2]+=(double)average(local_opt,num_local)/ratio;
			perf[3]+=(double)average(local_hemming,num_hemming)/ratio;
			
			if(loop==loop_num-Ts)
			{
				for(i=0;i<num_hemming;i++)fprintf(g,"%d %lf\n",d_hemming[i],local_hemming[i]);
				for(i=0;i<num_local;i++)fprintf(f,"%d %lf\n",d_local[i],local_opt[i]);
			}
		}
	}
	fclose(f);
	fclose(g);
	f=fopen("gnu-localopt","w");
	fprintf(f,"set term png\n");
	fprintf(f,"set output 'data/local-optimim-N:%d-K:%d-%d-%lf-%lf.png\n",N,K,rho,lambda,tau);
	fprintf(f,"plot '%s' pointsize 2 pointtype 6 title 'local optimum',",fn1);
	fprintf(f,"'%s' pointsize 2 pointtype 1 title 'local optimum with hemming distance'\n",fn2);
	fclose(f);
	system("gnuplot gnu-localopt");

	for(i=0;i<N;i++)
	{
		delete[] c[i];
		delete[] interaction[i];
	}
	delete[] c;
	delete[] interaction;
	delete[] global;
	delete[] global_future;
}

int org_forms(agent* a)
{
	int i,j;
	int *form;
	int num=0;
	form=new int[num_agent];
	for(i=0;i<num_agent;i++)
	{
		for(j=0;j<num;j++)
		{
			if(binary_2_number(a[i].seq)==form[j])break;
		}
		if(j==num)
		{
			form[num]=binary_2_number(a[i].seq);
			num++;
		}
	}
	delete[] form;
	return num;
}

int NK_timeseries(int rho,double lambda,double tau,int option)
{
        int i,j;
        double **c;
        int **interaction;
        agent *a;
        a=new agent[num_agent];
        for(i=0;i<num_agent;i++) a[i].sequence();

        c=new double*[N];
        for(i=0;i<N;i++)
	{
                c[i]=new double[n_k];
                for(j=0;j<n_k;j++) c[i][j]=uniform();
        }

        interaction=new int*[N];
        for(i=0;i<N;i++) interaction[i]=new int[K];
        interaction_allocation(interaction);

        int *global;
        global=new int[N];
	FILE *f;
	f=fopen("nk-time.csv","w");
	if(option==1) fprintf(f,"0 %d\n",org_forms(a));
	else if(option==2) fprintf(f,"0 %lf\n",a[0].cont);
        int* list;
        list=new int[N];
        for(i=0;i<N;i++) list[i]=i;
        for(int loop=1;loop<=loop_num;loop++)
	{
                if(loop%Ts==0)
		{
			reallocate_fitness(c,tau);
		}
                if(loop%Tg==0)
		{
                        get_global_optimum(global,c,interaction);
                        shuffle(list);
                        for(i=0;i<num_agent;i++)
			{
                                for(j=0;j<rho;j++)
				{
                                        if(uniform()<=lambda) a[i].seq[list[j]]=global[list[j]];
                                        else a[i].seq[list[j]]=1-global[list[j]];
                                }
                        }
                }
		for(i=0;i<num_agent;i++) a[i].local_search(c,interaction,list,rho);
                for(i=0;i<num_agent;i++) a[i].calculate_contribution(c,interaction);
		if(option==1) fprintf(f,"%d,%d\n",loop,org_forms(a));
		else if(option==2) fprintf(f,"%d %lf\n",loop,a[0].cont);
        }
	fclose(f);
	if(option==0)return org_forms(a);
	else if(option==2)
	{
		f=fopen("gnu-time","w");
		fprintf(f,"set title 'time series for firm 0 with N=%d,K=%d,rho=%d,lambda=%lf,tau=%lf'\n",N,K,rho,lambda,tau);
		fprintf(f,"plot 'nk-time.csv' with lines\n");
		fclose(f);
		system("gnuplot -persist gnu-time");
	}
        for(i=0;i<N;i++)
	{
                delete[] c[i];
                delete[] interaction[i];
        }
        delete[] c;
        delete[] interaction;
        delete[] a;
        delete[] global;
        delete[] list;
	return 0;
}

void NK_firmtimeseries(int rho,double lambda,double tau,int option)
{
        int i,j;
        double **c;
        int **interaction;
        agent *a;
        a=new agent[num_agent];
        for(i=0;i<num_agent;i++) a[i].sequence();

        c=new double*[N];
        for(i=0;i<N;i++)
	{
                c[i]=new double[n_k];
                for(j=0;j<n_k;j++) c[i][j]=uniform();
        }

        interaction=new int*[N];
        for(i=0;i<N;i++) interaction[i]=new int[K];
        interaction_allocation(interaction);
        int *global;
        global=new int[N];
	FILE *f;
	f=fopen("nk-firmtime.csv","w");
        int* list;
        list=new int[N];
	double* data;
	data=new double[loop_num];

        for(i=0;i<N;i++)list[i]=i;
        for(int loop=1;loop<=loop_num;loop++)
	{
		for(i=0;i<num_agent;i++) a[i].local_random(c,interaction);
                if(loop%Ts==0)
		{
			if(tau>0)
			{
				reallocate_fitness(c,tau);
			}
		}
                if(loop%Tg==0&&lambda>0)
		{
                        get_global_optimum(global,c,interaction);
                        shuffle(list);
                        for(i=0;i<num_agent;i++)
			{
                                for(j=0;j<rho;j++)
				{
                                        if(uniform()<=lambda) a[i].seq[list[j]]=global[list[j]];
                                        else a[i].seq[list[j]]=1-global[list[j]];
                                }
                        }
                }
                for(i=0;i<num_agent;i++)a[i].calculate_contribution(c,interaction);
		fprintf(f,"%d,%lf\n",loop,a[0].cont);
		data[loop]=a[0].cont;
        }
	fclose(f);
	printf("autocorrelation=%lf\n",autocorrelation(data,1,loop_num));
        for(i=0;i<N;i++)
	{
                delete[] c[i];
                delete[] interaction[i];
        }
        delete[] c;
        delete[] interaction;
        delete[] a;
        delete[] global;
        delete[] list;
}



int NK_sd_fitness(int rho,double lambda,double tau,int option)
{
	int i,j;
	double **c;
	int **interaction;
	agent *a;
	a=new agent[num_agent];
	for(i=0;i<num_agent;i++)a[i].sequence();

        c=new double*[N];
        for(i=0;i<N;i++)
	{
                c[i]=new double[n_k];
                for(j=0;j<n_k;j++) c[i][j]=uniform();
        }

        interaction=new int*[N];
        for(i=0;i<N;i++) interaction[i]=new int[K];
        interaction_allocation(interaction);

        int *global;
        global=new int[N];
	FILE *f;
	get_global_optimum(global,c,interaction);

	if(option==1)
	{
	        for(int loop=1;loop<=loop_num;loop++)
		{
			for(i=0;i<num_agent;i++)a[i].local_search2(c,interaction);
	                for(i=0;i<num_agent;i++)a[i].calculate_contribution(c,interaction);
	        }
		f=fopen("nk-sd.csv","w");
		for(i=0;i<num_agent;i++)
		{
			fprintf(f,"%d %lf\n",distance(a[i].seq,global),a[i].cont);
		}
		fclose(f);
	}
	else if(option==2)
	{
		int *temp;
		int size=(int)pow(2,N);
		temp=new int[N];
		f=fopen("nk-sd.csv","w");
		for(i=0;i<size;i++)
		{
			get_binary_from_number(i,temp);
			if(local_optimum(temp,c,interaction)) fprintf(f,"%d %lf\n",distance(global,temp),calculate_cont(temp,c,interaction));
		}
		fclose(f);
		f=fopen("gnu-sd","w");
		fprintf(f,"plot 'nk-sd.csv' pointsize 2 pointtype 6\n");
		fclose(f);
		system("gnuplot -persist gnu-sd");
	}
	else if(option==3||option==4)
	{
		int* list;
		list=new int[N];
		for(i=0;i<N;i++) list[i]=i;
		for(int loop=0;loop<loop_num;loop++)
		{
			if(loop%Ts==0)
			{
				reallocate_fitness(c,tau);
			}
			if(loop%Tg==0)
			{
				get_global_optimum(global,c,interaction);
				shuffle(list);
				for(i=0;i<num_agent;i++)
				{
					for(j=0;j<rho;j++)
					{
						if(uniform()<=lambda) a[i].seq[list[j]]=global[list[j]];
						else a[i].seq[list[j]]=1-global[list[j]];
					}
				}
			}

			for(i=0;i<num_agent;i++)a[i].local_search(c,interaction,list,rho);
			for(i=0;i<num_agent;i++)a[i].calculate_contribution(c,interaction);
		}
		int same=0;
		for(i=0;i<num_agent;i++)
		{
			for(j=0;j<N;j++)
			{
				if(global[j]!=a[i].seq[j]) break;
			}
			if(j==N) same++;
		}
		if(option==3)
		{
			f=fopen("nk-sd.csv","w");
			for(i=0;i<num_agent;i++)
			{
				fprintf(f,"%d %lf\n",distance(a[i].seq,global),a[i].cont);
			}
			fclose(f);
			f=fopen("gnu-sd","w");
			fprintf(f,"set title '# of agents that reach global optimum=%d' font \",15\"\n",same);
			fprintf(f,"set xlabel 'Strategic distance from global optimum' font \",15\"\n");
			fprintf(f,"set ylabel 'Firm value' font \",15\"\n");
			fprintf(f,"plot 'nk-sd.csv' with points lt -1 pointsize 2 pointtype 6 notitle\n");
			fclose(f);
			system("gnuplot -persist gnu-sd");
			char name[200];
			sprintf(name,"cp nk-sd.csv figure/sd-N_%d-K_%d-rho_%d-lambda_%.2lf-tau_%.2lf.csv",N,K,rho,lambda,tau);
			system(name);
			
			
		}
		else if(option==4)return same;
	}
        for(i=0;i<N;i++)
	{
                delete[] c[i];
                delete[] interaction[i];
        }
        delete[] c;
        delete[] interaction;
        delete[] a;
        delete[] global;
	return 0;
}
