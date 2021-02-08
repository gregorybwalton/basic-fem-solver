#include "fem.h"

double volume(double *a11,double *a21,double *a31,double *a41)
// calculate element volume using pointers for the each node point
// then increments for each component
{
        double vol;

        vol = 1.0*det3(*a21,*(a21+1),*(a21+2),*a31,*(a31+1),*(a31+2),*a41,*(a41+1),*(a41+2))\
                -*a11*det3(1.0,*(a21+1),*(a21+2),1.0,*(a31+1),*(a31+2),1.0,*(a41+1),*(a41+2))\
                +*(a11+1)*det3(1.0,*a21,*(a21+2),1.0,*a31,*(a31+2),1.0,*a41,*(a41+2))\
                -*(a11+2)*det3(1.0,*a21,*(a21+1),1.0,*a31,*(a31+1),1.0,*a41,*(a41+1));
        vol = vol/6.0;

        return vol;
}

double evolume(int el)
// just give it a element number
// easier then having to write it out each time
{
	double *s = msh.s;
	int *icon = msh.icon;
	int knode = msh.knode;
	int dim = msh.dim;
	double evol;

	evol = volume(&s[icon[el*knode]*dim],&s[icon[el*knode+1]*dim],&s[icon[el*knode+2]*dim],&s[icon[el*knode+3]*dim]);
	return evol;
}
