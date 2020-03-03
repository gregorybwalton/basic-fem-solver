#include "../fem.h"

void spsolve(spstiff,load,sol)
struct spmat spstiff;
double *load, *sol;
{
	int *bdflag = msh.bdflag;
	int nt = msh.ntrue;
	int nzeros = spstiff.nzeros;
	list_node **list = spstiff.head;
	list_node *current;
        int nnode = msh.nnode;
	int node, nbd;
	struct timeval start, end;
	struct rusage usage;
	double timet;

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
			current = list[nbd];
			while (current != NULL)
			{
				//printf("i = %d, j = %d, stiff = %.5f\n",nbd,current->na,current->val);
				gsl_spmatrix_set(K,nbd,current->na,current->val);
				current = current->next;
			}
			gsl_vector_set(f,nbd,load[nbd]);
			gsl_vector_set(u,nbd,sol[node]);
		}
	}

	// solve with GMRES iterative solver
	C = gsl_spmatrix_ccs(K);
	const double tol = 1.0e-6; // solution relative tolerance
	const size_t max_iter = 1000; // max iterations
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
			sol[node] = gsl_vector_get(u,nbd);
			//printf("sol(%.5f,%.5f,%.5f) = %.15f\n",msh.s[node][0],msh.s[node][1],msh.s[node][2],sol[node]);
		}
	}

	timet = (double)(end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) * 1e-6;
	printf("Time taken to solve (CPU time): %.5f s\n", timet);
        printf("--------------------\n");


}
