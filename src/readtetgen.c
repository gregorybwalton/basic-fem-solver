#include "fem.h"

void chkinterior(double *,int *);
void readnode();
void readele();
void readface();

mesh msh;
char ifnm[50];
char idir[50];
int soltype;

void readmesh(char* fnm,char* fdir)
// Reading in the meshes generated by TetGen
// fnm is the filename
// fdir is the file directory
// Needs be files with the extension .node and .ele
{
	int node,nintr;
	int abt;
	int nreg,n;

	strcpy(ifnm,fnm);
	strcpy(idir,fdir);

	printf("\n--------------------\n");
	printf("Reading in mesh...\n");

        readnode();
	readele();
	readface();

	nintr = 0; // Number of non-boundary nodes
	for (node=0;node<msh.nnode;node++)
	{
        	if (msh.bdflag[node]==0)
		{
			// interior
                	msh.bdflag[node] = nintr++; // should be post increment
	        }
	        else if (msh.bdflag[node]==1)
		{
			// boundary
        	        msh.bdflag[node] = -1;
	        }
		//printf("(%.15f,%.15f,%.15f), %i\n",msh.s[node*msh.dim+0],msh.s[node*msh.dim+1],msh.s[node*msh.dim+2],msh.bdflag[node]);
	}
	msh.neq = nintr;

	//msh.vol = (double *)malloc(sizeof(double)*msh.nel);

	printf("Number of dimensions = %d\n",msh.dim);
	printf("Number of nodes per element = %d\n",msh.knode);
	printf("Number of nodes = %d\n",msh.nnode);
	printf("Number of elements = %d\n",msh.nel);
	printf("Number of unknowns = %d\n",msh.neq);
	printf("--------------------\n");

}

void readnode( )
// Read .node file
{
	FILE *nfil;
	char nfnm[100];
	int nnode,dim,node,n,abt;
	int srt;

	snprintf(nfnm,sizeof(nfnm),"%s%s%s",idir,ifnm,".node");

	printf("Reading: %s\n",nfnm);

	if ((nfil = fopen(nfnm,"r"))==NULL)
	{
		printf("cannot open: %s\n",nfnm);
		exit(1);
	}

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
        	chkinterior(&msh.s[node*dim],&msh.bdflag[node]);
		//printf("(%.15f,%.15f,%.15f), %i\n",msh.s[node*dim+0],msh.s[node*dim+1],msh.s[node*dim+2],msh.bdflag[node]);
	}
        fclose(nfil);

	return;
}

void readele( )
// Read .ele file
{
	FILE *efil;
	char efnm[100];
	int nele,knode,nreg;
	int el,n;
	int srt;

	snprintf(efnm,sizeof(efnm),"%s%s%s",idir,ifnm,".ele");

	printf("Reading: %s\n",efnm);

	if ((efil = fopen(efnm,"r"))==NULL)
	{
		printf("cannot open: %s\n",efnm);
		exit(1);
	}

	srt = fscanf(efil,"%d %d %d\n",&nele,&knode,&nreg);

	msh.nel = nele;
	msh.knode = knode;
	msh.icon = (int *)malloc(sizeof(int)*nele*knode);
	if (nreg == 0)
	{
		msh.region = (int *)NULL;
	}
	else if (nreg == 1)
	{
		msh.region = (int *)malloc(sizeof(int)*nele);
	}

	for (el=0;el<nele;el++)
	{
		srt = fscanf(efil,"%*d");
		for (n=0;n<knode;n++)
		{
			srt = fscanf(efil," %d",&msh.icon[el*knode+n]);
			//printf("%d, ",msh.icon[el*knode+n]);
			//
			// Reduce index by one
			msh.icon[el*knode+n]--;
		}

		if (nreg == 1)
		{
			srt = fscanf(efil," %d",&msh.region[el]);
		}
		//printf("reg[%d] = %d\n",el,msh.region[el]);
		srt = fscanf(efil,"\n");
		//msh.region[el] = 0;

	}
	fclose(efil);

	return;
}

void readface(void)
// Read .face file
{
	FILE *ffil;
	char ffnm[100];
	int srt;
	int n,i,nf,bif,nface;
	int kf = msh.knode-1;
	
	snprintf(ffnm,sizeof(ffnm),"%s%s%s",idir,ifnm,".face");

	printf("Reading: %s\n",ffnm);

	if ((ffil = fopen(ffnm,"r"))==NULL)
	{
		printf("cannot open: %s\n",ffnm);
		exit(1);
	}
	
	srt = fscanf(ffil,"%d %d",&nface,&bif);
	//printf("nface = %d\n",nface);
	msh.nface = nface;
	
	msh.iconf = (int *)malloc(sizeof(int)*kf*nface);
	msh.bdflagf = (int *)malloc(sizeof(int)*nface);
	//msh.lscon = (double *)malloc(sizeof(double)*msh.nel);

	for (n=0;n<nface;n++)
	{
		srt = fscanf(ffil,"%*d %d %d %d",&msh.iconf[(n*kf)],&msh.iconf[(n*kf)+1],&msh.iconf[(n*kf)+2]);
		if (bif==1)
		{
			srt = fscanf(ffil," %d\n",&msh.bdflagf[n]);
			//printf("bdflagf[%d] = %d\n",n,(*bdflagf)[n]);
		}
		// reducing the index
		for (i=0;i<kf;i++)
		{
			msh.iconf[(n*kf)+i]--;
		}
	}
	fclose(ffil);
	return;
}

void chkinterior(double* s,int* bd)
{
	// check all interior nodes
	int ii;
	double tol = 1e-12;
	int chk = 1;
	int bdx1 = 1.0; // set the boundaries
	int bdx2 = 0.0;

	for (ii=0;ii<3;ii++)
	{
		if (fabs(*(s+ii)-bdx1)<tol || fabs(*(s+ii)-bdx2)<tol)
		{
			*bd = 1;
			//printf("(%.15f,%.15f,%.15f), %i\n",*(s+0),*(s+1),*(s+2),*bd);
			return;
		}
	}
	*bd = 0;
	//printf("(%.15f,%.15f,%.15f), %i\n",*(s+0),*(s+1),*(s+2),*bd);
	return;
}

void freemesh(mesh *tmsh)
{
	free(tmsh->icon);free(tmsh->s);free(tmsh->bdflag);free(tmsh->region);free(tmsh->lscon);

	return;
}
