#include "fem.h"

void chkinterior(double *,int *);

void readmesh(fnm,fdir)
char fnm[50];
char fdir[50];
// Reading in the meshes generated by TetGen
// fnm is the filename
// fdir is the file directory
// Needs be files with the extension .node and .ele
{
	char nfnm[100];
	char efnm[100];
	FILE *nfil;
	FILE *efil;
	int nnode,dim,nele,knode;
	int node,el,i,nintr;
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
        	chkinterior(&msh.s[node*dim],&msh.bdflag[node]);
		//printf("(%.15f,%.15f,%.15f), %i\n",msh.s[node*dim+0],msh.s[node*dim+1],msh.s[node*dim+2],msh.bdflag[node]);
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
			//printf("%d, ",msh.icon[el*knode+n]);
		}

		if (nreg == 1)
		{
			srt = fscanf(efil," %d",&msh.region[el]);
		}
		else
		{
			msh.region[el] = 0;
		}
		//printf("reg[%d] = %d\n",el,msh.region[el]);
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

	nintr = 0; // Number of non-boundary nodes
	for (node=0;node<nnode;node++)
	{
        	if (msh.bdflag[node]==0)
		{
			// interior
                	msh.bdflag[node] = nintr;
			nintr++;
	        }
	        else if (msh.bdflag[node]==1)
		{
			// boundary
        	        msh.bdflag[node] = -1;
	        }
		//printf("(%.15f,%.15f,%.15f), %i\n",msh.s[node*dim+0],msh.s[node*dim+1],msh.s[node*dim+2],msh.bdflag[node]);
	}
	msh.neq = nintr;

	printf("Number of dimensions = %d\n",msh.dim);
	printf("Number of nodes per element = %d\n",msh.knode);
	printf("Number of nodes = %d\n",msh.nnode);
	printf("Number of elements = %d\n",msh.nel);
	printf("Number of unknowns = %d\n",msh.neq);
	printf("--------------------\n");

}


void readvmesh(fnm,fdir)
// VMESH files from GRUMMP
// Not finished
char fnm[50];
char fdir[50];
{
        char nfnm[100];
        FILE *fil;
	int node;
        int srt;
        int nnode,nfaces,nbfaces,nverts,dim;

        snprintf(nfnm,sizeof(nfnm),"%s%s%s",fdir,fnm,".vmesh");
        printf("%s\n",nfnm);

        fil = fopen(nfnm,"r");
        srt = fscanf(fil,"%f %f %f %f\n",&nnode,&nfaces,&nbfaces,&nverts);
	dim = 3;

        for (node=0;node<nnode;node++)
        {
                srt = fscanf(fil,"%lf %lf %lf\n",&msh.s[node*dim+0],&msh.s[node*dim+1],&msh.s[node*dim+2]);
        }

}

void readvtk(fnm,fdir)
// VTK files from GRUMMP
char fnm[50];
char fdir[50];
{
        char nfnm[100];
	char buf[100];
        FILE *fil;
	int i,n,node,el,nei;
	int nnode,nele,knode,nint,dim;
	int verttyp;
	int srt;
	double *lensc;
	int *ctyp;
	
	printf("\n--------------------\n");
	printf("Reading in mesh...\n");

	strcpy(ifnm,fnm);
        snprintf(nfnm,sizeof(nfnm),"%s%s%s",fdir,fnm,".vtk");
        printf("%s\n",nfnm);

        if ((fil = fopen(nfnm,"r"))==NULL)
	{
		printf("cannot open %s\n",fnm);
		exit(1);
	}
	
	// header
        for (i=0;i<4;i++)
        {
                srt = fscanf(fil,"%*[^\n]\n");
                //fgets(buf,100,fil);
        }

	// nodes
        srt = fscanf(fil,"%*s %d %*s\n",&nnode);
	msh.nnode = nnode;
	dim = 3;
	msh.dim = dim;
	msh.s = (double *)malloc(sizeof(double)*nnode*dim);
	msh.bdflag = (int *)malloc(sizeof(int)*nnode);
        lensc = (double *)malloc(sizeof(double)*nnode);
	ctyp = (int *)malloc(sizeof(int)*nnode);
	for (node=0;node<nnode;node++)
        {
                srt = fscanf(fil,"%lf %lf %lf\n",&msh.s[node*dim+0],&msh.s[node*dim+1],&msh.s[node*dim+2]);
        }

	// elements
	srt = fscanf(fil,"%*s %d %*f\n",&nele);
	msh.nel = nele;
	nei = 0;
	for (el=0;el<nele;el++)
	{
		srt = fscanf(fil,"%d",&knode);
		if (el==0)
		{
			msh.knode = knode;
			msh.icon = (int *)malloc(sizeof(int)*nele*knode);
		}
		if (knode>3)
		{
			for (n=0;n<knode;n++)
			{
				srt = fscanf(fil," %d",&msh.icon[el*knode+n]);
			}
			nei++;
		}
		else
		{
			for (n=0;n<knode;n++)
			{
				srt = fscanf(fil," %*d");
			}
		}
		srt = fscanf(fil,"\n");
	}
	msh.nel = nei;
	//printf("msh.nel*knode = %d\n",msh.nel*knode);
	//msh.icon = (int *)realloc(msh.icon,sizeof(int)*msh.nel*knode);
	//exit(0);
	
	// cell types
	srt = fscanf(fil,"%*[^\n]\n");
	for (el=0;el<nele;el++)
	{
		srt = fscanf(fil,"%*f\n");
	}
	
	// vert type - not sure
        for (i=0;i<4;i++)
        {
                srt = fscanf(fil,"%*[^\n]\n");
        }
        
	for (node=0;node<nnode;node++)
        {
		srt = fscanf(fil,"%d\n",&verttyp);
		// Going to need to rely on checker
		chkinterior(&msh.s[node*dim],&msh.bdflag[node]);
		
		/*
		if (verttyp == 6)
		{
			msh.bdflag[node] = 0;
			nint++;
		}
		else
		{
			msh.bdflag[node] = 1;
		}
		*/
	}

	// scale length - useful for later
        for (i=0;i<2;i++)
        {
                srt = fscanf(fil,"%*[^\n]\n");
        }
        for (node=0;node<nnode;node++)
        {
		srt = fscanf(fil,"%lf\n",&lensc[node]);
	}

	// Cell type - is this regions?
        for (i=0;i<3;i++)
        {
                srt = fscanf(fil,"%*[^\n]\n");
        }
	msh.region = (int *)malloc(sizeof(int)*msh.nel);
        for (el=0;el<msh.nel;el++)
        {
		//srt = fscanf(fil,"%lf\n",&ctyp[node]);
		// watch out this could be a problem - inconsistent with the previous length of region
		srt = fscanf(fil,"%d\n",&msh.region[el]);
		msh.region[el]--;
		//printf("region[%d] = %d\n",el,msh.region[el]);
	}
	fclose(fil);

        int nintr = 0; // Number of non-boundary nodes
        for (node=0;node<nnode;node++)
        {
                if (msh.bdflag[node]==0)
                {
                        // interior
			msh.bdflag[node] = nintr;
                        nintr++;
		}
                else if (msh.bdflag[node]==1)
                {
                	// boundary
                        msh.bdflag[node] = -1;
                }
	}
        msh.neq = nintr;

        printf("Number of dimensions = %d\n",msh.dim);
        printf("Number of nodes per element = %d\n",msh.knode);
        printf("Number of nodes = %d\n",msh.nnode);
        printf("Number of elements = %d\n",msh.nel);
        printf("Number of unknowns = %d\n",msh.neq);
        printf("--------------------\n");
}

void chkinterior(s,bd)
double *s;
int *bd;
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
