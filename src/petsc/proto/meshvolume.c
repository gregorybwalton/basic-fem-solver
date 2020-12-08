#include "../petsc.h"

void writevols(double *);
void halfrefinement(double *,int *);

void meshvolume( )
{
	double *volmax = (double *)calloc(msh.nel,sizeof(double));
	double *vol = (double *)calloc(msh.nel,sizeof(double));

	printf("RUNNING MESH REFINEMENT\n");

	calcvol(vol);
	halfrefinement(volmax,vol);

	writevols(volmax);
	return;
}

void calcvol(vol)
double *vol;
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

void internalrefinement(volmax,vol)
double *volmax,*vol;
{
	int el,i,j;
	int knode = msh.knode;
	int dim = msh.dim;
	double *cent = (double *)calloc(dim,sizeof(double));

	for (el=0;el<msh.nel;el++)
	{
		for (i=0;i<dim;dim++)
		{
			cent[i] = 0.0;
			for (j=0;j<knode;j++)
			{
				cent[i] += msh.s[msh.icon[el*knode+j]*dim+i];
			}
			cent[i] /= knode;
		}
		

		// need to indentify the surface mesh
		for (i=0;i<dim;dim++)
		{
		}		


	}
}

void halfrefinement(volmax,vol)
double *volmax,*vol;
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
			volmax[el] = kappa*vol[el];
		}
		else
		{
			volmax[el] = vol[el]; // neg means uncontrained
		}
	}
	return;
}

void writevols(volmax)
double *volmax;
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
