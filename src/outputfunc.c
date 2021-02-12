#include "fem.h"

char sfnm[50];

char *savename(void)
{
	char tfnm[50];
	char *sptr = sfnm;
	if (soltype==0)
	{
		strcpy(tfnm,"linear");
	}
	else if (soltype==1)
	{
		strcpy(tfnm,"quad");
	}
	//sprintf(sfnm,"%s/%s_%d","./results",tfnm,msh.nel); // problem type and number of elements
	sprintf(sfnm,"%s/%s","./results",ifnm);
	
	return sptr;
}

void resoutput(double* res)
{
        int nnode = msh.nnode;
        //double **s = msh.s;
        double *s = msh.s;
	int dim = msh.dim;
        int node;
	char rfnm[50];
        FILE *resfil;
        
        strcpy(rfnm,savename());
        strcat(rfnm,".res");

	printf("Saving: %s\n",rfnm);
        resfil = fopen(rfnm,"w");
	fprintf(resfil,"%d\n",nnode);
        for (node=0;node<nnode;node++)
        {
                //fprintf(resfil,"%.15f %.15f %.15f %.15f\n",s[i][0],s[i][1],s[i][2],res[i]);
                fprintf(resfil,"%.15f %.15f %.15f %.15f\n",s[node*dim+0],s[node*dim+1],s[node*dim+2],res[node]);
        }
        fclose(resfil);
}

void erroroutput(double nerror,double* error)
{
        int nnode = msh.nnode;
        double *s = msh.s;
        int dim = msh.dim;
        int node;
        FILE *efil;
        char efnm[50];

        strcpy(efnm,savename());
        strcat(efnm,".err");

	printf("Saving: %s\n",efnm);
        efil = fopen(efnm,"w");
	fprintf(efil,"%d\n",nnode);
        fprintf(efil,"%.15f\n",nerror);
        for (node=0;node<nnode;node++)
        {
                fprintf(efil,"%.15f %.15f %.15f %.15f\n",s[node*dim+0],s[node*dim+1],s[node*dim+2],error[node]);
        }
        fclose(efil);
}
      

void writevol()
{
	FILE *fil;
	int el;
	int nel = msh.nel;
	int *reg = msh.region;

	double vol;
	double tvol = 0.;
	double invol = 0.;
	for (el=0;el<nel;el++)
	{
		vol = evolume(el);
		tvol += vol;
		if (reg[el]==1) invol += vol;
	}

	fil = fopen("results/vols.out","a");
	fprintf(fil,"%.15f, %.15f\n",tvol,invol);

	fclose(fil);

	return;
}
