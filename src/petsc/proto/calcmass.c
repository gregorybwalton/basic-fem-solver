#include "../../fem.h"

double evol(int);

double calcmass(void)
// calcuting the Me*u for each element and summing
// total mesh mass
{
	int nel = msh.nel;
	int knode = msh.knode;
	int *icon = msh.icon;
	int el,i,j,jj;
	double vol,emass,tmass;
	
	tmass = 0.0;
	for (el=0;el<nel;el++)
	{
		vol = evol(el)*(1.0/20.0);
		emass = 0.0;
		for (i=0;i<knode;i++)
		{
			for (j=0;j<knode;j++)
			{
				jj = icon[el*knode+j];
				if (i==j)
				{
					emass += 2.0*vol*spstiff.sol[jj];
				}
				else
				{
					emass += vol*spstiff.sol[jj];
				}
			}
		}
		tmass += emass;
	}
	
	return tmass;
}

double evol(int el)
// calculate element volume - using same formulate in d3n4.c
{
	int i,j;
	double vol;
	double *s = msh.s;
	int *icon = msh.icon;
	int knode = msh.knode;
	int dim = msh.dim;
	double *sloc = (double *) calloc(knode*dim,sizeof(double));

        for (i=0;i<knode;i++)
	{
                for (j=0;j<dim;j++)
		{
        		sloc[i*dim+j] = s[icon[el*knode+i]*dim+j];
		}
	}
	
	vol = det3(sloc[3],sloc[4],sloc[5],sloc[6],sloc[7],sloc[8],sloc[9],sloc[10],sloc[11])\
				-sloc[0]*det3(1.0,sloc[4],sloc[5],1.0,sloc[7],sloc[8],1.0,sloc[10],sloc[11])\
                                +sloc[1]*det3(1.0,sloc[3],sloc[5],1.0,sloc[6],sloc[8],1.0,sloc[9],sloc[11])\
                                -sloc[2]*det3(1.0,sloc[3],sloc[4],1.0,sloc[6],sloc[7],1.0,sloc[9],sloc[10]);
	vol = (1.0/6.0)*fabs(vol);
	
	free(sloc);
	return vol;
}
