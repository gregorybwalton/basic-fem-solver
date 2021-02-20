#include "../../fem.h"

void halfrefinevol(double *);
void internalrefinevol();
void writevols(double *);

// useful for when need to specify the volume refinement requirement at each cell
void meshvolumescale( )
{
	// do I really need vol and lscon, can I just one and rewrite
	double *vol = (double *)calloc(msh.nel,sizeof(double));
	msh.lscon = (double *)calloc(msh.nel,sizeof(double));

	printf("Running mesh refinement\n");

	//halfrefinevol(vol);
	internalrefinevol();

	//writevols(msh.lscon);
	return;
}

void internalrefinevol()
// Refinement given the distance to an interior boundary
{
	int el,i,j,n;
	int knode = msh.knode;
	int dim = msh.dim;
	int nel = msh.nel;
	double *s = msh.s;
	int *icon = msh.icon;
	int *iconf = msh.iconf;
	int *bdff = msh.bdflagf;
	int nface = msh.nface;
	double vol;
	double r,p;
	
	double *cent = (double *)calloc(dim,sizeof(double));
	double *rmin = (double *)calloc(nel,sizeof(double));

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
		rmin[el] = 1000.;
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
					if ((r<rmin[el]) && (r!=0.0)) rmin[el] = sqrt(r);
					//if ((r<rmin[el]) && (r!=0.0)) rmin[el] = r;
				}	
			}
		}
		//printf("rmin[%d] = %.10e\n",el,rmin[el]);
		vol = volume(&s[icon[el*knode]*dim],&s[icon[el*knode+1]*dim],&s[icon[el*knode+2]*dim],&s[icon[el*knode+3]*dim]);
		//printf("vol = %.5e\n",vol);
		msh.lscon[el] = rmin[el]*vol;
		//printf("lscon[%d] = %.5e\n",el,msh.lscon[el]);
	}

	free(cent); free(rmin);
	return;
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
