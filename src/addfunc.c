#include "fem.h"

struct gauss preGauss(void)
{
	struct gauss gaussval;
	int knode = msh.knode;
	//printf("preGauss knode = %d\n",knode);
        gaussval.w = (double *)calloc(knode,sizeof(double));
        gaussval.crd = (double *)calloc(knode*knode,sizeof(double));
        double alpha, beta;

	//printf("knode = %i\n",knode);
	if (knode == 4)
	{
        alpha = 0.58541020;
        beta = 0.13819660;
        gaussval.w[0] = 1.0/4.0;
        gaussval.w[1] = gaussval.w[0];
        gaussval.w[2] = gaussval.w[0];
        gaussval.w[3] = gaussval.w[0];
        gaussval.crd[0] = alpha;
        gaussval.crd[1] = beta;
        gaussval.crd[2] = beta;
        gaussval.crd[3] = beta;
        gaussval.crd[4] = beta;
        gaussval.crd[5] = alpha;
        gaussval.crd[6] = beta;
        gaussval.crd[7] = beta;
        gaussval.crd[8] = beta;
        gaussval.crd[9] = beta;
        gaussval.crd[10] = alpha;
        gaussval.crd[11] = beta;
        gaussval.crd[12] = beta;
        gaussval.crd[13] = beta;
        gaussval.crd[14] = beta;
        gaussval.crd[15] = alpha;
	}
	else if (knode == 3)
	{
        gaussval.w[0] = 1.0/3.0;
        gaussval.w[1] = 1.0/3.0;
        gaussval.w[2] = 1.0/3.0;
        gaussval.crd[0] = 1.0/2.0;
        gaussval.crd[1] = 1.0/2.0;
        gaussval.crd[2] = 0.0;
        gaussval.crd[3] = 0.0;
        gaussval.crd[4] = 1.0/2.0;
        gaussval.crd[5] = 1.0/2.0;
        gaussval.crd[6] = 1.0/2.0;
        gaussval.crd[7] = 0.0;
        gaussval.crd[8] = 1.0/2.0;
	}

	return gaussval;
}


double calcGauss(s,node,gaussvals)
double *s;
int node;
struct gauss gaussvals;
{
        int i;
        double detJ, csp;
        double *w;
	double *coord;
	int knode = msh.knode;
	
	w = gaussvals.w;
	coord = gaussvals.crd;
	
	detJ = det3(s[3]-s[0],s[6]-s[0],s[9]-s[0],s[4]-s[1],s[7]-s[1],\
			s[10]-s[1],s[5]-s[2],s[8]-s[2],s[11]-s[2]);
        csp = 1.0/6.0; // Volume of unit tetrahedron

        double xzeta, yzeta, zzeta;
        double f; // loading function
        double Nbar = 0.0;
        double F = 0.0;
        for (i=0;i<knode;i++)
        {
                xzeta = s[0] + (s[3]-s[0])*coord[0*knode+i]+(s[6]-s[0])*coord[1*knode+i]+(s[9]-s[0])*coord[2*knode+i];
                yzeta = s[1] + (s[4]-s[1])*coord[0*knode+i]+(s[7]-s[1])*coord[1*knode+i]+(s[10]-s[1])*coord[2*knode+i];
                zzeta = s[2] + (s[5]-s[2])*coord[0*knode+i]+(s[8]-s[2])*coord[1*knode+i]+(s[11]-s[2])*coord[2*knode+i];
                f = loadFunc(xzeta,yzeta,zzeta);
		//printf("loadFunc = %.15f\n",f);

                switch(node)
                {
                        case 0:
                                Nbar = 1.0 - coord[0*knode+i] - coord[1*knode+i] - coord[2*knode+i];
                        case 1:
                                Nbar = coord[0*knode+i];
                        case 2:
                                Nbar = coord[1*knode+i];
                        case 3:
                                Nbar = coord[2*knode+i];
                }


                F = F + (w[i]*f*Nbar*detJ*csp);
		/*
		printf("w[%i] = %.15f\n",i,w[i]);
		printf("f = %.15f\n",f);
		printf("detJ = %.15f\n",detJ);
		printf("Nbar = %.15f\n",Nbar);
		printf("csp = %.15f\n",csp);
                printf("F = %.15f\n",F);
                exit(0);
		*/
        }

        return F;
}

double det3(a11,a12,a13,a21,a22,a23,a31,a32,a33)
double a11,a12,a13,a21,a22,a23,a31,a32,a33;
{
// Calculate the determinate of 3 x 3 matrix
        double det = a11*(a22*a33-a32*a23)\
                -a12*(a21*a33-a31*a23)\
                +a13*(a21*a32-a31*a22);

        return det;
}

double phiFunc(x,y,z)
double x,y,z;
{
// Define solution function for the boundary
        double phi;
        if (soltype==0)
        {
                phi = x + y + z; // linear
        }
        else if (soltype==1)
        {
                phi = x*(x-1.0)*y*(y-1.0)*z*(z-1.0); // quad
        }
	else
	{
		printf("Error: provide valid solution type");
		phi = 0.0;
	}

        return phi;
}

double loadFunc(x,y,z)
double x,y,z;
{
// Define the loading function
        double f;
        if (soltype==0)
        {
                f = 0.0;
        }
        else if (soltype==1)
        {
                f = -2.0*(y*(y-1.0)*z*(z-1.0)+x*(x-1.0)*z*(z-1.0)+x*(x-1.0)*y*(y-1.0));
        }
	else
	{
		printf("Error: provide valid solution type");
		f = 0.0;
	}
        return f;
}

