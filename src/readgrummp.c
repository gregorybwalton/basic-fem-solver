#include "fem.h"

void chkinterior(double *,int *);

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
	int verttyp,knmin;
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
	dim = 3; // Will always be 3d, but for 2d, x3=0.0
	knmin = 3; // The minimum number of nodes per element to be read
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
	nei = 0;
	for (el=0;el<nele;el++)
	{
		srt = fscanf(fil,"%d",&knode);
		if (el==0)
		{
			msh.knode = knode;
			icontmp = (int *)malloc(sizeof(int)*nele*knode);
		}
		//printf("knode[%d] = %d\n",el,knode);
		if (knode>knmin)
		{
			nei++;
		}
		for (n=0;n<knode;n++)
		{
			srt = fscanf(fil," %d",&icontmp[el*knode+n]);
		}
		srt = fscanf(fil,"\n");
	}
	msh.nel = nei;
	nels = nele;
	msh.icon = (int *)malloc(sizeof(int)*msh.nel*knode);
	for (el=0;el<msh.nel;el++)
	{
		for (n=0;n<knode;n++)
		{
			msh.icon[el*knode+n]=icontmp[el*knode+n];
		}
	}
	
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
