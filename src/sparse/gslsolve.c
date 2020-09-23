#include "../fem.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_splinalg.h>

void writeres(char *,int,int,int,int,double,double,double);
double calcerror(double *,int);

void spsolve(spstiff)
struct sysmat spstiff;
{
	int *bdflag = msh.bdflag;
	int nt = msh.neq;
	int nzeros = spstiff.nzeros;
        int nnode = msh.nnode;
	int nel = msh.nel;
	int dim = msh.dim;
	list_node *current;
	int node, nbd;
	struct timeval start, end;
	struct rusage usage;
	double timet,nerr;

	// In form Ku=f
	gsl_spmatrix *C;
	gsl_spmatrix *K = gsl_spmatrix_alloc_nzmax(nt,nt,nzeros,GSL_SPMATRIX_COO);
	gsl_vector *f = gsl_vector_alloc(nt);
	gsl_vector *u = gsl_vector_alloc(nt);

        printf("Solving sparse iteratively...\n");

	// populate gsl sparse matrix format
	for (node=0;node<nnode;node++)
	{
		nbd = bdflag[node];
		if (nbd>=0)
		{
			current = spstiff.head[nbd];
			//printf("row %d:",nbd);
			while (current != NULL)
			{
				gsl_spmatrix_set(K,nbd,current->na,current->val);
				//printf(" (%d, %g)",current->na,current->val);	

				current = current->next;
			}
			//printf("\n");
			gsl_vector_set(f,nbd,spstiff.load[nbd]);
			gsl_vector_set(u,nbd,spstiff.sol[node]);
			//printf("b[%i] = %g\n",nbd,spstiff.load[nbd]);
			//printf("x[%i] = %g\n",nbd,spstiff.sol[node]);
		}
	}
	

	// solve with GMRES iterative solver
	C = gsl_spmatrix_ccs(K);
	const double tol = 1.0e-6; // solution relative tolerance
	const size_t max_iter = 10000; // max iterations
	const gsl_splinalg_itersolve_type *T = gsl_splinalg_itersolve_gmres;
	gsl_splinalg_itersolve *work = gsl_splinalg_itersolve_alloc(T, nt, 0);
	size_t iter = 0;
	double residual;
	int status;

	getrusage(RUSAGE_SELF, &usage);
	start = usage.ru_utime;
	printf("  iter  |     residual  \n");
	do
	{
		status = gsl_splinalg_itersolve_iterate(C, f, tol, u, work);
		residual = gsl_splinalg_itersolve_normr(work);
		printf("%5zu %23.12e\n", iter, residual);

		if (status == GSL_SUCCESS) printf("Converged.\n");
	}
	while (status == GSL_CONTINUE && ++iter < max_iter);

	if (iter == max_iter) printf("WARNING: Did not converge!\n");

	getrusage(RUSAGE_SELF, &usage);
	end = usage.ru_utime;

	// get solution
        for (node=0;node<nnode;node++)
        {
                nbd = bdflag[node];
                if (nbd>=0)
		{
			spstiff.sol[node] = gsl_vector_get(u,nbd);
			printf("sol(%.5f,%.5f,%.5f) = %.15f\n",msh.s[node*dim+0],msh.s[node*dim+1],msh.s[node*dim+2],spstiff.sol[node]);
		}
	}

	timet = (double)(end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) * 1e-6;
	nerr = calcerror(spstiff.sol,1);
	writeres("results/rdata",nt,nnode,nel,iter,timet,residual,nerr);
	printf("Time taken to solve (CPU time): %.5f s\n", timet);
        printf("--------------------\n");
}

void writeres(wfil,ukn,nn,nel,its,timet,norm,err)
char wfil[50];
int ukn,nn,nel,its;
double timet,norm,err;
{
        FILE *fp;
        if (access(wfil,W_OK) != -1)
        {
                fp = fopen(wfil,"a");
                fprintf(fp,"%i,%i,%i,%i,%g,%g,%g\n",ukn,nn,nel,its,timet,norm,err);
                fclose(fp);
        }
}

