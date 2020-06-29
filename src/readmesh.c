#include "fem.h"

void readmesh(fnm,fdir)
char fnm[50];
char fdir[50];
{
	char nfnm[100];
	char efnm[100];
	FILE *nfil;
	FILE *efil;
	int nnode,dim;
	int node,el,i;
	int nele,knode;
	int abt;
	int nreg,n;
	int srt;

	strcpy(ifnm,fnm);

	printf("\n--------------------\n");
	printf("Reading in mesh...\n");

	snprintf(nfnm,sizeof(nfnm),"%s%s%s",fdir,ifnm,".node");

	printf("%s\n",nfnm);

	nfil = fopen(nfnm,"r");

	srt = fscanf(nfil,"%d %d %*d %*d\n",&nnode,&dim);
	
	msh.nnode = nnode;
	msh.dim = dim;
	msh.s = (double *)malloc(sizeof(double)*nnode*dim);
	msh.bdflag = (int *)malloc(sizeof(int)*nnode);

	for (node=0;node<nnode;node++)
	{
		srt = fscanf(nfil,"%*d");
		for (n=0;n<dim;n++)
		{
			srt = fscanf(nfil," %lf",&msh.s[node*dim+n]);
		}
		srt = fscanf(nfil," %d\n",&abt);
	                
		msh.bdflag[node] = abs(abt);
        	//printf("s(:,%d) = %.2f %.2f %.2f, bdflag = %d\n",node,msh.s[node*dim],msh.s[node*dim+1],msh.s[node*dim+2],msh.bdflag[node]);
	}
        fclose(nfil);
        
	snprintf(efnm,sizeof(efnm),"%s%s%s",fdir,ifnm,".ele");

	printf("%s\n",efnm);

	efil = fopen(efnm,"r");
	srt = fscanf(efil,"%d %d %d\n",&nele,&knode,&nreg);
	msh.nel = nele;
	msh.knode = knode;
	msh.icon = (int *)malloc(sizeof(int)*nele*knode);
	msh.region = (int *)malloc(sizeof(int)*nele);

	for (el=0;el<nele;el++)
	{
		srt = fscanf(efil,"%*d");
		for (n=0;n<knode;n++)
		{
			srt = fscanf(efil," %d",&msh.icon[el*knode+n]);
		}

		if (nreg == 1)
		{
			srt = fscanf(efil," %d",&msh.region[el]);
		}
		else
		{
			msh.region[el] = 0;
		}
		//printf("region[%d] = %d\n",el,msh.region[el]);
		srt = fscanf(efil,"\n");
		//msh.region[el] = 0;

		// Reduce the index by one
		//printf("icon[%d] = ",el);	
		for (i=0;i<knode;i++)
		{
			msh.icon[el*knode+i]--;
			//printf("%d, ",msh.icon[el*knode+i]);
		}
		//printf("\n");
	}
	fclose(efil);

	int nintr = 0; // Number of non-boundary nodes
	for (node=0;node<nnode;node++)
	{
        	if (msh.bdflag[node]==0)
		{
                	msh.bdflag[node] = nintr;
			nintr++;
	        }
	        else if (msh.bdflag[node]==1)
		{
        	        msh.bdflag[node] = -1;
	        }
		//printf("msh bdflag[%d] = %d\n",node,msh.bdflag[node]);
	}
	msh.neq = nintr;

	printf("Number of dimensions = %d\n",msh.dim);
	printf("Number of nodes per element = %d\n",msh.knode);
	printf("Number of nodes = %d\n",msh.nnode);
	printf("Number of elements = %d\n",msh.nel);
	printf("Number of unknowns = %d\n",msh.neq);
	printf("--------------------\n");

}
