#include "../fem.h"

void printspmatrices(spstiff,load)
struct sysmat spstiff;
double *load;
{
	int nt = msh.neq;
	
	printf("Writing stiffness matrix and load vector to file.\n");
        FILE *kfil;
        kfil = fopen("k_matrix","w");
        FILE *lfil;
        lfil = fopen("f_vector","w");

	list_node *current;
        int i;
        for (i=0;i<nt;i++)
        {
                fprintf(lfil,"%d %.10f\n",i,load[i]);
		current = spstiff.head[i];
		while (current != NULL)
		{
			fprintf(kfil,"%d %d %.10f\n",i,current->na,current->val);
			current = current->next;
                }
        }
        fclose(kfil);
        fclose(lfil);
}
