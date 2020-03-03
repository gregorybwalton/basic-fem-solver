#include "fem.h"

void errorouput(double,double *);

void calcerror(u)
double *u;
{
	printf("Calculating the error...\n");
	int nnode = msh.nnode;
	double *s = msh.s;
	int dim = msh.dim;
	int node;
	double uex;
	double *error;
	double nerror;

	error = (double *)calloc(nnode,sizeof(double));
	for (node=0;node<nnode;node++)
	{
		uex = phiFunc(s[node*dim+0],s[node*dim+1],s[node*dim+2]);
		error[node] = fabs(u[node]-uex);
		nerror = nerror + error[node]*error[node];
	}
	nerror = sqrtf(nerror)/nnode;
	printf("||Error|| = %.15f\n",nerror);

	erroroutput(nerror,error);

	free(error);
}

