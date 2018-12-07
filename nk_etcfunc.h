int get_number_from_binary(int x, int *seq, int** inter)
{
	int i,z=1,ret=0;
	for(i=K-1;i>=0;i--)
	{
		if(seq[inter[x][i]]==1) ret+=z;
		z*=2;
	}
	if(seq[x]==1) ret+=z;
	return ret;
}

void get_binary_from_number(int x, int* seq)
{
        int t=x,i;
        for(i=0;i<N;i++)
	{
                seq[N-i-1]=t%2;
                t/=2;
        }
}

double uniform()
{
        return (double)rand()/RAND_MAX;
}

int redundancy(int *seq,int x,int l)
{
        int i,j;
        for(i=0;i<l;i++)
	{
                for(j=i+1;j<l;j++)
		{
                        if(seq[i]==seq[j]) return 1;
                }
        }
        if(seq[l-1]==x) return 1;
        return 0;
}

void shuffle(int* list)
{
        int x,y;
        int temp;
        for(int loop=0;loop<N*20;loop++)
	{
                x=rand()%N;
                y=rand()%N;
                temp=list[x];
                list[x]=list[y];
                list[y]=temp;
        }
}

int binary_2_number(int* a)
{
        int i,j=1;
        int num=0;
        for(i=N-1;i>=0;i--)
	{
                num+=a[i]*j;
                j*=2;
        }
        return num;
}

double autocorrelation(double* data, int timelag, int N)
{
	int i;
	double *x1,*x2;
	double avr1=0,avr2=0,var1=0,var2=0,autocor=0;
	x1=new double[N-timelag];
	x2=new double[N-timelag];
	for(i=0;i<N-timelag;i++)
	{
		x1[i]=data[i];
		x2[i]=data[i+timelag];
	}
	for(i=0;i<N-timelag;i++)
	{
		avr1+=x1[i]/(N-timelag);
		avr2+=x2[i]/(N-timelag);
		var1+=x1[i]*x1[i]/(N-timelag);
		var2+=x2[i]*x2[i]/(N-timelag);
	}
	var1-=avr1*avr1;
	var2-=avr2*avr2;
	for(i=0;i<N-timelag;i++) autocor+=(x1[i]-avr1)*(x2[i]-avr2)/(N-timelag);
	return autocor/sqrt(var1*var2);
}

int fixed_decision(int x,int *fix,int rho)
{
	int i;
	for(i=0;i<rho;i++)
	{
		if(x==fix[i])return 1;
	}
	return 0;
}

void print(int* a)
{
	for(int i=0;i<N;i++) printf("%d ",a[i]);
}

double hemming_distance(int* a, int* b, int N)
{
	int i;
	double hemming=0;
	for(i=0;i<N;i++)
	{
		if(a[i]==b[i]) hemming++;
	}
	return hemming;
}

double average(double *a, int x)
{
	double avr=0;
	for(int i=0;i<x;i++) avr+=a[i]/x;
	return avr;
}
