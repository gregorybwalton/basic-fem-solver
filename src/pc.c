#include "fem.h"

void jacobi(spstiff,sol)
// Jacobi preconditioner, simply divide through by the leading diagonal
struct spmat spstiff;
double *sol;
{
	printf("Jacobi preconditioner.\n");
	int nnode = msh.ntrue;
	list_node **list = spstiff.head;
	list_node *current;
	int node,nbd;
	double diag;

	for (node=0;node<nnode;node++)
	{
		current = list[node];

		diag = current->val;
		while (current != NULL)
		{
			current->val = current->val/diag;
			current = current->next;
		}
		load[node] = load[node]/diag;
	}

	//printlist(list,nnode,1);
}

void ichol(spstiff,sol)
// Incomplete Cholesky factorisation
// Follows the algorithm as presented in "INCOMPLETE CHOLESKY FACTORIZATIONS WITH LIMITED
// MEMORY" by Lin and More
//
// WARNING: From testing, it only produces a correct lower diagonal matrix
struct spmat spstiff;
double *sol;
{
	printf("Incomplete Cholesky factorisiation.\n");

        int nn = msh.ntrue;
        int i, j, k;
	list_node **list = spstiff.head;
	int nnz = spstiff.nzeros;	
	list_node *aij,*aik,*ajk;

	writelist(list,nnz,nn,"list");
	
        for (j=0;j<nn;j++)
        {
		list[j]->val = sqrt(list[j]->val); // square the diagonals
		/*
		if (list[j]->val < 0)
		{
			list[j]->val = -sqrt(fabs(list[j]->val));
		}
		else
		{
			list[j]->val = sqrt(list[j]->val);
		}
		*/

		// Will need to sort the linked lists to improve efficiency
		// Subtracting the strict lower diagonal by column product
		for (k=0;k<j;k++)
		{
			ajk = list[j];
	
			while (ajk != NULL)
			{
				if (ajk->na == k) break;
				ajk = ajk->next;
			}

			for (i=j+1;i<nn;i++)
			{
				aij = list[i];
				aik = list[i];
				while (aij != NULL)
				{
					if (aij->na == j) break;
					aij = aij->next;
				}

				while (aik != NULL)
				{
					if (aik->na == k) break;
					aik = aik->next;
				}
				
				if (ajk != NULL & aij != NULL & aik != NULL)
				{
					printf("i = %d, j = %d, k = %d\n",i,j,k);
					//printf("aij = %.5f, aik = %.5f, ajk = %.5f\n",aij->val,aik->val,ajk->val);
					aij->val -= aik->val*ajk->val;
					//printf("aij = %.5f\n",aij->val);
				}
			}
		}

		// Dividing by Lij, and subtracting the diagonal by Lij^2
		for (i=j+1;i<nn;i++)
		{
			aij = list[i];
			while (aij != NULL)
			{
				if (aij->na == j) break;
				aij = aij->next; 
			}
			
			if (aij != NULL)
			{
				aij->val = aij->val/list[j]->val;
				list[i]->val -= aij->val*aij->val;
			}
		}
	}

	printf("Removing upper diagonal...\n");
	list_node *curr;
	// remove the upper diagonal, since its all junk
	for (j=0;j<nn;j++)
	{
		for (curr = list[j]; curr != NULL; curr = curr->next)
		{
			if (curr->na > j)
			{
				listdel(&list[j],curr->na);
			}
		}
	}

        writelist(list,nnz,nn,"listicf");
	// Now multiply the LHS and RHS by the factored matrix

	printf("BREAKING!\n");
	exit(0);
}
