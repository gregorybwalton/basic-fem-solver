#include "../../fem.h"
//#include "tetgen.h"

void inputmsh(tetgenio &);
void outputmsh(tetgenio &);
void interpsol(tetgenio &);
double volume(double *,double *,double *,double *);

void runtetgen()
{
	printf("\nTETGEN: runtetgen\n");
	tetgenio in,out;
	//char switches[50] = "rkqaO0"; // switches being used in refinement
	tetgenbehavior switches;
	switches.refine = 1; // -r, refines
	switches.vtkview = 0; // -k, outputs vtks, this doesn't seem to do anything
	switches.quality = 1; // -q, quality
	switches.varvolume = 1; // -a, varvolume
	switches.opt_iterations = 0;
	switches.opt_scheme = 0;
	switches.opt_max_flip_level = 0;
	switches.quiet =  1; // supresses terminal output, except errors

	in.firstnumber = 0;
	inputmsh(in);
	
	// addin, -i, additional vertices
	// bgmin, -m, constraint functions
	tetrahedralize(&switches,&in,&out);
	
	//out.save_nodes("example_refine");
	interpsol(out);
	outputmsh(out);

	void deinitialize();
	return;
}

void inputmsh(tetgenio &in)
{
	int n,d,el,k;
	int dim = msh.dim;
	int knode = msh.knode;
	int nel = msh.nel;

	in.mesh_dim = dim;
	in.numberofpoints = msh.nnode;
	//in.pointlist = new REAL[msh.nnode * msh.dim];
	in.pointlist = (REAL *)malloc(sizeof(REAL)*msh.nnode*msh.dim);
	in.numberofpointattributes = 1;
	//in.pointattributelist = new REAL[msh.nnode];
	in.pointattributelist = (REAL *)malloc(sizeof(REAL)*msh.nnode);
	
        for (n=0;n<msh.nnode;n++)
        {
                for (d=0;d<dim;d++)
                {
			in.pointlist[dim*n+d] = (REAL)msh.s[dim*n+d];
		}
		if (msh.bdflag[n] > -1)
                {
			in.pointattributelist[n] = 0.0;
		}
		else
		{
			in.pointattributelist[n] = -1.0;
		}

	}

	in.numberoftetrahedra = msh.nel;
	//in.tetrahedronlist = new int[msh.nel*msh.knode];
	in.tetrahedronlist = (int *)malloc(sizeof(int)*msh.nel*msh.knode);
	//in.tetrahedronvolumelist = new REAL[msh.nel];
	in.tetrahedronvolumelist = (REAL *)malloc(sizeof(REAL)*msh.nel);
	in.numberoftetrahedronattributes = 0;
	if (msh.region != NULL)
	{
		in.numberoftetrahedronattributes = 1;
		//in.tetrahedronattributelist = new REAL[msh.nel];
		in.tetrahedronattributelist = (REAL *)malloc(sizeof(REAL)*msh.nel);
	}

        for (el=0;el<nel;el++)
        {
		in.tetrahedronvolumelist[el] = (REAL)msh.lscon[el];
                for (k=0;k<knode;k++)
                {
                        in.tetrahedronlist[knode*el+k] = msh.icon[knode*el+k];
                }
		if (in.numberoftetrahedronattributes == 1) in.tetrahedronattributelist[el] = (REAL)msh.region[el];
        }

	return;
}

// Apparently this usage of pointer in the function declaration is unique to C++
void outputmsh(tetgenio &out)
{
	int n,d,el,k;
	int knode = msh.knode;
	msh.nnode = (int)out.numberofpoints;
	msh.nel = (int)out.numberoftetrahedra;

	printf("TETGEN: Outputting new mesh\n");

	msh.s = (double*)realloc(msh.s,sizeof(double)*msh.dim*msh.nnode);
	msh.bdflag = (int*)realloc(msh.bdflag,sizeof(int)*msh.nnode);
	int nintr = 0;
	for (n=0;n<msh.nnode;n++)
	{
		for (d=0;d<msh.dim;d++)
		{
			msh.s[msh.dim*n+d] = (double)out.pointlist[msh.dim*n+d];
		}
		if (out.pointattributelist[n]==0) 
		{
			msh.bdflag[n] = nintr++;
		}
		else if (out.pointattributelist[n]!=0)  
		{
			msh.bdflag[n] = -1;
		}
	}

	msh.icon = (int*)realloc(msh.icon,sizeof(int)*msh.nel*msh.knode);
	if (out.numberoftetrahedronattributes == 1) msh.region = (int*)realloc(msh.region,sizeof(int)*msh.nel);
	for (el=0;el<msh.nel;el++)
	{
		for (k=0;k<msh.knode;k++)
		{
			msh.icon[knode*el+k] = (int)out.tetrahedronlist[knode*el+k];
		}
		if (out.numberoftetrahedronattributes == 1) msh.region[el] = (int)out.tetrahedronattributelist[el];
		//printf("el = %d\n",el);
	}

	return;
}

void interpsol(tetgenio &out)
{
	int n,el,i,j,inck;
	double v,vk;
	int knode = msh.knode;
	int dim = msh.dim;
	double *newu = (double *)malloc(sizeof(double)*out.numberofpoints);
	double *beta = (double *)malloc(sizeof(double)*knode);
	double **ik = (double **)malloc(sizeof(double *)*knode);
	int count = 0; // using for debugging

	printf("TETGEN: Interpolating to new mesh\n");

	for (n=0;n<out.numberofpoints;n++) newu[n] = 1.0;

	for (n=0;n<out.numberofpoints;n++)
	{
		for (el=0;el<msh.nel;el++)
		{
			/*
			printf("p =\n");
			printf("[%.5f %.5f %.5f;\n",msh.s[msh.icon[knode*el]*dim],msh.s[msh.icon[knode*el]*dim+1],msh.s[msh.icon[knode*el]*dim+2]);
			printf("%.5f %.5f %.5f;\n",msh.s[msh.icon[knode*el+1]*dim],msh.s[msh.icon[knode*el+1]*dim+1],msh.s[msh.icon[knode*el+1]*dim+2]);
			printf("%.5f %.5f %.5f;\n",msh.s[msh.icon[knode*el+2]*dim],msh.s[msh.icon[knode*el+2]*dim+1],msh.s[msh.icon[knode*el+2]*dim+2]);
			printf("%.5f %.5f %.5f];\n",msh.s[msh.icon[knode*el+3]*dim],msh.s[msh.icon[knode*el+3]*dim+1],msh.s[msh.icon[knode*el+3]*dim+2]);
			*/
			v = volume(&msh.s[msh.icon[knode*el]*dim],&msh.s[msh.icon[knode*el+1]*dim],&msh.s[msh.icon[knode*el+2]*dim],&msh.s[msh.icon[knode*el+3]*dim]);
			//printf("vol[%d] = %.5e\n",el,v);
			inck = 1;
			for (i=0;i<knode;i++)
			{
				for (j=0;j<knode;j++) 
				{
					ik[j] = &(msh.s[msh.icon[knode*el+j]*dim]);
				}
				ik[i] = &(out.pointlist[n*dim]);
				//printf("pk[%d] = [%.5f %.5f %.5f];\n",i,*ik[i],*(ik[i]+1),*(ik[i]+2));
				vk = volume(ik[0],ik[1],ik[2],ik[3]);
				beta[i] = vk/v;
				// dodgy rounding - round to 10 d.p. since the vtk is to 6 d.p
				double tmp = (long int)(beta[i]*1e10 + .5);
				beta[i] = tmp / 1e10;

				//printf("beta[%d] = %.5f\n",i,beta[i]);
				//if (fabs(beta[i]) < 1.e-5 && beta[i]!=0.0) printf("beta[%d] = %.5e\n",i,beta[i]);
				/*
				if ((fabs(out.pointlist[n*dim]-0.65)<1.e-6 && fabs(out.pointlist[n*dim+1]-0.35)<1.e-6 && fabs(out.pointlist[n*dim+2]-0.4625)<1.e-6) && i==(knode-1))
				if ((fabs(out.pointlist[n*dim]-0.65)<1.e-6 && fabs(out.pointlist[n*dim+1]-0.545878)<1.e-6 && fabs(out.pointlist[n*dim+2]-0.610504)<1.e-6 && fabs(msh.s[msh.icon[knode*el]*dim]-0.65)<1.e-6))
				{
					//printf("[%.5f %.5f %.5f]\n",out.pointlist[n*dim],out.pointlist[n*dim+1],out.pointlist[n*dim+2]);
					printf("	beta[%d] = %.5e\n",i,beta[i]);
					printf("pk[%d] = [%.5f %.5f %.5f];\n",i,*ik[i],*(ik[i]+1),*(ik[i]+2));
					printf("p =\n");
					printf("[%.5f %.5f %.5f;\n",msh.s[msh.icon[knode*el]*dim],msh.s[msh.icon[knode*el]*dim+1],msh.s[msh.icon[knode*el]*dim+2]);
					printf("%.5f %.5f %.5f;\n",msh.s[msh.icon[knode*el+1]*dim],msh.s[msh.icon[knode*el+1]*dim+1],msh.s[msh.icon[knode*el+1]*dim+2]);
					printf("%.5f %.5f %.5f;\n",msh.s[msh.icon[knode*el+2]*dim],msh.s[msh.icon[knode*el+2]*dim+1],msh.s[msh.icon[knode*el+2]*dim+2]);
					printf("%.5f %.5f %.5f];\n",msh.s[msh.icon[knode*el+3]*dim],msh.s[msh.icon[knode*el+3]*dim+1],msh.s[msh.icon[knode*el+3]*dim+2]);
				}
				*/
				//if (fabs(beta[i]) < 1.e-10 && beta[i] != 0.0) printf("beta[%d] = %.5e\n",i,beta[i]);

				if (beta[i] < 0.0)
				{
					inck = 0;
					break;
				}

			}
			if (inck == 1)
			{
				newu[n] = 0.0;
				// interpolation based upon volume ratio
				for (i=0;i<knode;i++)
				{
					newu[n] += beta[i]*spstiff.sol[msh.icon[knode*el+i]];
				}
				count++;
				//if (count==5) printf("	newu[%d] = %.5f\n",n,newu[n]);
			}
		}
	}
	spstiff.sol = (double*)realloc(spstiff.sol,sizeof(double)*out.numberofpoints);
	for (n=0;n<out.numberofpoints;n++)
	{
		spstiff.sol[n] = newu[n];
	}
	free(newu);

	return;
}

double volume(double *a11,double *a21,double *a31,double *a41)
{
	double vol;

	vol = 1.0*det3(*a21,*(a21+1),*(a21+2),*a31,*(a31+1),*(a31+2),*a41,*(a41+1),*(a41+2))\
		-*a11*det3(1.0,*(a21+1),*(a21+2),1.0,*(a31+1),*(a31+2),1.0,*(a41+1),*(a41+2))\
		+*(a11+1)*det3(1.0,*a21,*(a21+2),1.0,*a31,*(a31+2),1.0,*a41,*(a41+2))\
		-*(a11+2)*det3(1.0,*a21,*(a21+1),1.0,*a31,*(a31+1),1.0,*a41,*(a41+1));
	vol = vol/6.0;

	return vol;
}
