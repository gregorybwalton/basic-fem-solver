#include "../petsc.h"


void writels(double *);
void internalrefinels(double *,int *);
void halfrefinels(double *,int *);

void meshlengthscale( )
{
	//int *icon = icontmp;
	//internalrefinement(rmin,icon);

	halfrefinels(msh.lscon,msh.icon);

	writels(msh.lscon);
	return;
}

void internalrefinels( )
{
	int i,j,m,n;
	int iel,jel;
	int knode = msh.knode;
	int nel = msh.nel;
	double r;
	
	double *rmin = (double *)calloc(msh.nnode,sizeof(double));
	
	for (iel=0;iel<nel;iel++)
	{
		for (n=0;n<knode;n++)
		{
			rmin[msh.icon[iel*knode+n]] = 100.0;
			//printf("mesh.nel=%d, nels=%d\n",msh.nel,nels);
			for (jel=nel;jel<nels;jel++)
			{
				for (m=0;m<3;m++)
				{
					r = 0.0;
					for (i=0;i<3;i++)
					{
						r += pow(msh.s[msh.icon[iel*knode+n]+i]-msh.s[msh.icon[jel*3+m]+i],2);
					}
					r = pow(r,0.5);
					if ((r < rmin[msh.icon[iel*knode+n]]) && r != 0.0)
					{
						rmin[msh.icon[iel*knode+n]] = r*r;
					}
				}
			}
			//printf("r[%d] = %.5f\n",iel*4+n,rmin[iel*4+n]);
			msh.lscon[msh.icon[iel*knode+n]] *= rmin[msh.icon[iel*knode+n]];
		}
	}
	return;
}

void halfrefinels( )
{
	int el,n;
	int knode = msh.knode;
	int nel = msh.nel;

	for (el=0;el<nel;el++)
	{
		for (n=0;n<knode;n++)
		{
			if (msh.s[msh.icon[el*knode+n]]<0.5)
			{
				msh.lscon[msh.icon[el*knode+n]] *= 1.0e-1;
			}
			//else
			//{
			//	rmin[icon[el*knode+n]] = 1.e-1;
			//}
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
