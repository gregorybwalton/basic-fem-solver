#include "../fem.h"

double easpectratio(int);
double calclength(double *,double *,int);
double calcheight(double *,double *,double *,double *);

void meshquality(double *minvol,double *maxaspectratio)
{
	int nel = msh.nel;
	int el;
	double ar,vol;

	double maxar = -1.e6;
	double minar = -1.*maxar;
	double maxv = -100.;
	double minv = -1.*maxv;

	for (el=0;el<nel;el++)
	{
		ar = easpectratio(el);
		vol = evolume(el);
		if (ar>maxar) maxar = ar;
		if (ar<minar) minar = ar;
		if (vol>maxv) maxv = vol;
		if (vol<minv) minv = vol;

	}

	printf("Cell aspect ratio\n");
	printf("Max, min: %3e, %3e\n",maxar,minar);
	printf("Cell volume\n");
	printf("Max, min: %3e, %3e\n",maxv,minv);

	if (minv<0.0) printf("WARNING: negative volumes!\n");
	
	if (minvol) *minvol = minv;
	if (maxaspectratio) *maxaspectratio = maxar;

	return;
}


double easpectratio(int el)
// calculate the aspect ratio for each element
// max length over the minimum height
{
	int i,j;
	double *s = msh.s;
	int *icon = msh.icon;
	int dim = msh.dim;
	int knode = msh.knode;
	double elen, ear;
	double h;

	double lmax = -1.e3;
	double lmin = 1.e3;
	double hmin = 1.e3;
	for (i=0;i<(knode-1);i++)
	{
		for (j=i+1;j<knode;j++)
		{
			elen = calclength(&s[icon[el*knode+i]*dim],&s[icon[el*knode+j]*dim],dim);

			if (elen < lmin) lmin = elen;
			if (elen > lmax) lmax = elen;
		}
	}
	for (i=0;i<knode;i++)
	{
		h = calcheight(&s[icon[el*knode+(i%knode)]*dim],&s[icon[el*knode+((i+1)%knode)]*dim],&s[icon[el*knode+((i+2)%knode)]*dim],&s[icon[el*knode+((i+3)%knode)]*dim]);
		if (h<hmin) hmin = h;
	}

	//ear = lmax/lmin;
	ear = lmax/hmin; // tetgen aspect ratio, see manual 

	return ear;
}

double calclength(double *xi,double *xj,int dim)
{
	double len = 0.0;
	int i;

	for (i=0;i<dim;i++)
	{
		len += (*(xi+i)-*(xj+i))*(*(xi+i)-*(xj+i));
	}
	len = sqrt(len);

	return len;
}

double calcheight(double *p1, double *p2, double *p3, double *p4)
{
	int i;
	double d,vw,vm;
	double v[3],w[3];

	// cross product
	v[0] = ((*(p2+1)-*(p1+1))*(*(p3+2)-*(p1+2))) - ((*(p2+2)-*(p1+2))*(*(p3+1)-*(p1+1)));
	v[1] = ((*(p2+2)-*(p1+2))*(*(p3)-*(p1))) - ((*(p2)-*(p1))*(*(p3+2)-*(p1+2)));
	v[2] = ((*(p2)-*(p1))*(*(p3+1)-*(p1+1))) - ((*(p2+1)-*(p1+1))*(*(p3)-*(p1)));

	//printf("v = %.4f, %.4f, %.4f\n",v[0],v[1],v[2]); // correct
	vw = 0.0;
	vm = 0.0;
	for (i=0;i<3;i++) 
	{
		w[i] = *(p1+i)-*(p4+i);
		vw += v[i]*w[i];
		vm += v[i]*v[i];
		
	}
	vw = fabs(vw);
	vm = sqrt(vm);
	d = vw/vm;
	//printf("w = %.4f, %.4f, %.4f\n",w[0],w[1],w[2]); // correct
	//printf("vw = %.4f, vm = %.4f\n",vw,vm);
	return d;
}

// to do
void skewness(void)
{
	return;
}
