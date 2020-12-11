#include "petsc.h"

void calcvol(double *);
void halfrefinevol(double *);
void internalrefinevol(double *);
int readface(int **,int **);
void writevols(double *);

void meshvolume( )
{
	// do I really need vol and lscon, can I just one and rewrite
	double *vol = (double *)calloc(msh.nel,sizeof(double));
	msh.lscon = (double *)calloc(msh.nel,sizeof(double));

	printf("RUNNING MESH REFINEMENT\n");

	calcvol(vol);
	//halfrefinevol(vol);
	internalrefinevol(vol);

	writevols(msh.lscon);
	return;
}

void calcvol(double* vol)
{
	int el,n,i;
	int knode = msh.knode;
	int dim = msh.dim;

	int nv = knode-1; // number of vectors

	double *a = (double *)calloc(dim*nv,sizeof(double));

	for (el=0;el<msh.nel;el++)
        {
		for (n=0;n<nv;n++)
		{
			for (i=0;i<dim;i++)
			{
				a[n*dim+i] =  msh.s[msh.icon[el*knode+n+1]*dim+i]-msh.s[msh.icon[el*knode]*dim+i];
				//printf("a[%d] = %.5e\n",n*nv+i,a[n*nv+i]);
			}
		}
		vol[el] = (1.0/6.0)*fabs(a[0]*(a[4]*a[8]-a[5]*a[7])+a[1]*(a[5]*a[6]-a[3]*a[8])+a[2]*(a[3]*a[7]-a[4]*a[6])); // vol of each element
		//printf("vol[%d] = %.5e\n",el,vol[el]);
	}

	return;
}

void internalrefinevol(double* vol)
// Refinement given the distance to an interior boundary
// requires a .face file
{
	int el,i,j,n;
	int knode = msh.knode;
	int dim = msh.dim;
	int nel = msh.nel;
	double *cent;
	double r,*rmin;
	double p;
	

	cent = (double *)calloc(dim,sizeof(double));
	rmin = (double *)calloc(nel,sizeof(double));

	// read in the face data
	int nface;
	int *iconf,*bdff;
	nface = readface(&iconf,&bdff);

	// finding all the interior boundary nodes
	// could be used to make it more efficient
	/*
	int nfaceinner = 0;
	double *iconfi;
	for (n=0;n<nface;n++)
	{
		//printf("bdff[%d] = %d\n",n,bdff[n]);
		if (bdff[n]=0)
		{
			nfaceinner += 1;
			realloc(iconfi,nfaceinner*(knode-1));
			for (i=0;i<(knode-1);i++)
			{
				iconfi[nfaceinner*(knode-1)+i]=iconf[n*(knode-1)+i];
			}
		}
	}
	*/

	for (el=0;el<msh.nel;el++)
	{
		for (i=0;i<dim;i++)
		{
			cent[i] = 0.0;
			for (j=0;j<knode;j++)
			{
				//printf("msh.icon[%d] = %d\n",el*knode+j,msh.icon[el*knode+j]);
				cent[i] += msh.s[msh.icon[el*knode+j]*dim+i];
			}
			cent[i] /= knode;
		}
		rmin[el] = 100.0;	
		for (n=0;n<nface;n++)
		{
			//printf("bdff[%d] = %d\n",n,bdff[n]);
			if (bdff[n]==0)
			{
				//printf("bdff[%d] = %d\n",n,bdff[n]);
				for (i=0;i<(knode-1);i++)
				{
					r = 0.0;
					for (j=0;j<dim;j++)
					{
						p = msh.s[iconf[n*(knode-1)+i]*dim+j];
						//printf("%.5e, ",p);
						r += pow(cent[j] - p,2.0);
					}
					//printf("\n");
					//r = pow(r,0.5); // just for a square					
					if ((r<rmin[el]) && (r!=0.0)) rmin[el] = r;
				}	
			}
		}
		//printf("rmin[%d] = %.10e\n",el,rmin[el]);
		msh.lscon[el] = rmin[el]*vol[el];
	}
}

void halfrefinevol(double* vol)
{
	int el,n;
	int knode = msh.knode;
	int dim = msh.dim;

	double xcent; // x centroid
	double kappa = 1.0e-2;
	for (el=0;el<msh.nel;el++)
	{
		xcent = 0.0;
		for (n=0;n<knode;n++)
		{
			xcent += msh.s[msh.icon[el*knode+n]*dim];
		}
		xcent = xcent/knode;
		//printf("xcent[%d] = %.5f\n",el,xcent);
		if (xcent<0.5)
		{
			msh.lscon[el] = kappa*vol[el];
		}
		else
		{
			msh.lscon[el] = vol[el]; // neg means uncontrained
		}
	}
	return;
}

void writevols(double* volmax)
{
	FILE *fil;
	int el;
	char vfnm[50];

	snprintf(vfnm,sizeof(vfnm),"%s%s",ifnm,".vol");

	if ((fil = fopen(vfnm,"w")) == NULL)
	{
		perror("fopen()");
		exit(1);
	}
	fprintf(fil,"%d\n",msh.nel);

	for (el=0;el<msh.nel;el++)
	{
		// don't want to refine on the boundary
		//if (msh.bdflag[node] != -1)
		//{
			// if its non-zero
			if (volmax[el] != 0.0)
			{
				fprintf(fil,"%5d%22.16f\n",el+1,volmax[el]);
			}
		//}
	}
	fclose(fil);
	return;
}
