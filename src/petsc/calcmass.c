#include "../fem.h"

double calcmass(void)
// calcuting the Me*u for each element and summing
// total mesh mass
{
	int nel = msh.nel;
	int knode = msh.knode;
	int dim = msh.dim;
	double *s = msh.s;
	int *icon = msh.icon;
	int el,i,j,jj;
	double vol,emass,tmass;
	
	tmass = 0.0;
	for (el=0;el<nel;el++)
	{
		vol = volume(&s[icon[el*knode]*dim],&s[icon[el*knode+1]*dim],&s[icon[el*knode+2]*dim],&s[icon[el*knode+3]*dim])*(1./20.);
		emass = 0.0;
		for (i=0;i<knode;i++)
		{
			for (j=0;j<knode;j++)
			{
				jj = icon[el*knode+j];
				if (i==j)
				{
					emass += 2.0*vol*spstiff.sol[jj];
				}
				else
				{
					emass += vol*spstiff.sol[jj];
				}
			}
		}
		tmass += emass;
	}
	
	return tmass;
}

