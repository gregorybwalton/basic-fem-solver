#include "../petsc.h"

void writels(double *);
void internalrefinement(double *,int *);
void halfrefinement(double *,int *);

void meshlengthscale( )
{
	double *rmin = (double *)calloc(msh.nnode,sizeof(double));

	//int *icon = icontmp;
	//internalrefinement(rmin,icon);

	halfrefinement(rmin,msh.icon);

	writels(rmin);
	return;
}

void internalrefinement(double *rmin,int *icon)
{
	int i,j,m,n;
	int iel,jel;
	int knode = msh.knode;
	double r;
	
	for (iel=0;iel<msh.nel;iel++)
	{
		for (n=0;n<knode;n++)
		{
			if (rmin[icon[iel*knode+n]]==0.0)
			{
				rmin[icon[iel*knode+n]] = 100.0;
			}
			//printf("mesh.nel=%d, nels=%d\n",msh.nel,nels);
			for (jel=msh.nel;jel<nels;jel++)
			{
				for (m=0;m<3;m++)
				{
					r = 0.0;
					for (i=0;i<3;i++)
					{
						r += pow(msh.s[icon[iel*knode+n]+i]-msh.s[icon[jel*3+m]+i],2);
					}
					r = pow(r,0.5);
					if ((r < rmin[icon[iel*knode+n]]) && r != 0.0)
					{
						rmin[icon[iel*knode+n]] = r*r*r;
					}
				}
			}
			//printf("r[%d] = %.5f\n",iel*4+n,rmin[iel*4+n]);
		}
	}
	return;
}

void halfrefinement(double *rmin,int *icon)
{
	int el,n;
	int knode = msh.knode;

	for (el=0;el<msh.nel;el++)
	{
		for (n=0;n<knode;n++)
		{
			if (msh.s[icon[el*knode+n]]<0.5)
			{
				rmin[icon[el*knode+n]] = 1.e-2;
			}
			else
			{
				rmin[icon[el*knode+n]] = 1.e-1;
			}
		}
	}
	return;
}

void writels(double *rmin)
{
	FILE *fil;
	int node;

	fil = fopen("lenscale","w");

	for (node=0;node<msh.nnode;node++)
	{
		// don't want to refine on the boundary
		//if (msh.bdflag[node] != -1)
		//{
			// if its non-zero
			if (rmin[node] != 0.0)
			{
				fprintf(fil,"%d %4.10e\n",node,rmin[node]);
			}
		//}
	}
	fclose(fil);
	return;
}
