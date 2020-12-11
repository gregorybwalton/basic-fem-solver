#include "fem.h"

void erroroutput(double,double *);

double calcerror(double* u,int pidx)
{
	printf("Calculating the error...\n");
	int nnode = msh.nnode;
	double *s = msh.s;
	int dim = msh.dim;
	int node;
	double uex;
	double *abse; // absolute error
	double *rele; // relative error
	double nabse = 0.0;
	double nrele = 0.0;

	abse = (double *)calloc(nnode,sizeof(double));
	rele = (double *)calloc(nnode,sizeof(double));
	for (node=0;node<nnode;node++)
	{
		uex = phiFunc(s[node*dim+0],s[node*dim+1],s[node*dim+2]);
		abse[node] = fabs(u[node]-uex);
		rele[node] = abse[node]/uex;
		nabse = nabse + abse[node]*abse[node];
		nrele = nrele + rele[node]*rele[node];
	}
	//nabse = sqrt(nabse)/nnode;
	nabse = sqrt(nabse);
	printf("||Error|| = %.15f\n",nabse/nnode);

	if (pidx == 1)
	{
		erroroutput(nabse,abse);
	}

	free(abse);
	free(rele);
	return nabse;
}

