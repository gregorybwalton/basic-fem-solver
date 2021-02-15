#include "../fem.h"
#include "tetgen.h"

void runtetgen();
void meshvolumescale();
void inputmsh(tetgenio &);
void outputmsh(tetgenio &,mesh &);
void interpsol(double *,int);
void interplocal(double *,int);
double calcmass();
void tetfin(tetgenio &);

void meshquality(double *,double *);

void checkqualtetgen()
// check the mesh quality and determine a re mesh
{
	double minvol,maxar;

	meshquality(&minvol,&maxar);

	if (minvol < 1.e-8 || maxar > 1.e2) runtetgen();

	return;
}

void runtetgen()
// running re meshing
{
	printf("\nTETGEN: runtetgen\n");
	
	//void initialize();
	tetgenio in,out;
	in.firstnumber = 0; // first index is 0
	tetgenbehavior switches;

	// define switches
	switches.refine = 1; // -r, refines
	switches.quality = 1; // -q, quality
	switches.neighout = 1; // -n, neighout
	//switches.nobisect = 1; // -Y, supresses boundary facets/segments splitting
	//switches.diagnose = 1; // -d, detecting self intersections
	if (msh.lscon != NULL) switches.varvolume = 1; // -a, varvolume
	switches.opt_max_asp_ratio = 1.e2; // max aspect ratio
	//switches.opt_iterations = 1; // max number of iterations
	//switches.opt_scheme = 0; // local mesh operations, 0 to 7
	//switches.opt_max_flip_level = 0; // max flip level, 0 to 10
	//switches.minratio = 4.0; // -q, min ratio
	switches.quiet = 1; // supresses terminal output, except errors
	switches.verbose = 0; // 0 nothing, max 4

	double totalm1 = calcmass();

	inputmsh(in); // converts to tetgen mesh format
	// addin, -i, additional vertices
	// bgmin, -m, constraint functions
	tetrahedralize(&switches,&in,&out);
	//out.save_nodes("testout");
	//interpsol(out.pointlist,out.numberofpoints);
	interplocal(out.pointlist,out.numberofpoints);
	outputmsh(out,msh);
	//tetfin(in);
	void deinitialize(); // don't need tetgen anymore

        double totalm2 = calcmass();
        double percentdm = ((totalm1 - totalm2)/totalm1)*100.0;
        printf("Change in total mass = %.5f %\n",percentdm);

	return;
}

void tetfin(tetgenio &in)
{
	free(in.pointlist);
	free(in.tetrahedronlist);
	free(in.pointmarkerlist);
	if (msh.lscon != NULL) free(in.tetrahedronvolumelist);
	free(in.trifacelist);
	free(in.trifacemarkerlist);
	if (msh.region != NULL) free(in.tetrahedronattributelist);

	return;
}

void inputmsh(tetgenio &in)
{
	int n,d,el,k;
	int dim = msh.dim;
	int knode = msh.knode;
	int nel = msh.nel;
	int nnode = msh.nnode;
	int nface = msh.nface;

	in.mesh_dim = dim;
	in.numberoftetrahedra = nel;
	in.numberofpoints = nnode;
	in.numberoftrifaces = nface;
	in.numberofpointattributes = 0;
	in.numberoftetrahedronattributes = 0;

	// 'new' style C++ dynamic allocation
	in.pointlist = new REAL [nnode*dim];
	in.tetrahedronlist = new int [nel*knode];
	in.pointmarkerlist = new int [nnode];
	if (msh.lscon != NULL) in.tetrahedronvolumelist = new REAL [nel];
	in.trifacelist = new int [(knode-1)*nface];
	in.trifacemarkerlist = new int [nface];

	// C dynamic allocation
	/*
	in.pointlist = (REAL *)malloc(sizeof(REAL)*nnode*dim);
	in.tetrahedronlist = (int *)malloc(sizeof(int)*nel*knode);
	in.pointmarkerlist = (int *)malloc(sizeof(int)*nnode);
	if (msh.lscon != NULL) in.tetrahedronvolumelist = (REAL *)malloc(sizeof(REAL)*nel);
	in.trifacelist = (int *)malloc(sizeof(int)*(knode-1)*nface);
	in.trifacemarkerlist = (int *)malloc(sizeof(int)*nface);
	*/
	if (msh.region != NULL)
	{
		in.numberoftetrahedronattributes = 1;
		in.tetrahedronattributelist = new REAL[msh.nel];
		//in.tetrahedronattributelist = (REAL *)malloc(sizeof(REAL)*msh.nel);
	}
	
        for (n=0;n<nnode;n++)
        {
                for (d=0;d<dim;d++)
                {
			in.pointlist[dim*n+d] = (REAL)msh.s[dim*n+d];
		}
		if (msh.bdflag[n] > -1)
                {
			in.pointmarkerlist[n] = 0;
		}
		else
		{
			in.pointmarkerlist[n] = -1;
		}

	}

        for (el=0;el<nel;el++)
        {
		if (msh.lscon != NULL) in.tetrahedronvolumelist[el] = (REAL)msh.lscon[el];
                for (k=0;k<knode;k++)
                {
                        in.tetrahedronlist[knode*el+k] = msh.icon[knode*el+k];
                }
		if (in.numberoftetrahedronattributes == 1) in.tetrahedronattributelist[el] = (REAL)msh.region[el];
        }

	for (n=0;n<nface;n++)
	{
                for (k=0;k<(knode-1);k++)
                {
			in.trifacelist[(knode-1)*n+k] = msh.iconf[(knode-1)*n+k];
                }
		in.trifacemarkerlist[n] = msh.bdflagf[n];
        }

	return;
}

// Apparently this usage of pointer in the function declaration is unique to C++
void outputmsh(tetgenio &out,mesh &msho)
// the second argument can just be the original mesh(?)
{
	int n,i,el,k;
	int knode = msh.knode;
	int dim = msh.dim;
	int nnode = (int)out.numberofpoints;
	int nel = (int)out.numberoftetrahedra;
	int nface = (int)out.numberoftrifaces;

	// update new structure
	msho.nnode = nnode;
	msho.nel = nel;
	msho.dim = dim;
	msho.knode = knode;
	msho.nface = nface;

	printf("TETGEN: Outputting new mesh\n");
	
	if (&msho == &msh)
	{
		msho.icon = (int*)realloc(msho.icon,sizeof(int)*nel*knode);
		msho.s = (double*)realloc(msho.s,sizeof(double)*dim*nnode);
		msho.bdflag = (int*)realloc(msho.bdflag,sizeof(int)*nnode);
		if (out.numberoftetrahedronattributes == 1) msho.region = (int*)realloc(msho.region,sizeof(int)*nel);
		msho.iconf = (int*)realloc(msho.iconf,sizeof(int)*nface*(knode-1));
		msho.bdflagf = (int*)realloc(msho.bdflagf,sizeof(int)*nface);
		msho.neigh = (int *)realloc(msho.neigh,sizeof(int)*nel*knode);
	}
	else
	{
		msho.icon = (int*)malloc(sizeof(int)*nel*knode);
		msho.s = (double*)malloc(sizeof(double)*dim*nnode);
		msho.bdflag = (int*)malloc(sizeof(int)*nnode);
		if (out.numberoftetrahedronattributes == 1) msho.region = (int*)malloc(sizeof(int)*nel);
		msho.iconf = (int*)malloc(sizeof(int)*nface*(knode-1));
		msho.bdflagf = (int*)malloc(sizeof(int)*nface);
		msho.neigh = (int*)malloc(sizeof(int)*nel*knode);
	}

	int nintr = 0;
	double x,y,z;
	double ep=1.e-8;
	for (n=0;n<nnode;n++)
	{
		for (i=0;i<dim;i++)
		{
			msho.s[dim*n+i] = (double)out.pointlist[dim*n+i];
		}
		if (out.pointmarkerlist[n]==0) 
		{
			msho.bdflag[n] = nintr++;
		}
		else if (out.pointmarkerlist[n]==1 || out.pointmarkerlist[n]==-1)  
		{
			msho.bdflag[n] = -1;
			//printf("%d\n",out.pointmarkerlist[n]);
			
			// botch: tetgen seems to be return some interior nodes as -1
			// may need more investigation
			x = msho.s[dim*n];
			y = msho.s[dim*n+1];
			z = msho.s[dim*n+2];
			if (fabs(x)>ep && fabs(x-1.)>ep && fabs(y)>ep && fabs(y-1.)>ep && fabs(z)>ep && fabs(z-1.)>ep)
			{
				//printf("(%.5e,%.5e,%.5e), %d\n",x,y,z,out.pointmarkerlist[n]);
				msho.bdflag[n] = nintr++;
			}	
		}

	}
	msho.neq = nintr;
	
	for (el=0;el<nel;el++)
	{
		for (k=0;k<knode;k++)
		{
			msho.icon[knode*el+k] = (int)out.tetrahedronlist[knode*el+k];
			msho.neigh[knode*el+k] = (int)out.neighborlist[knode*el+k];
		}
		if (out.numberoftetrahedronattributes == 1) msho.region[el] = (int)out.tetrahedronattributelist[el];
			//printf("el = %d\n",el);
	}

	for (n=0;n<nface;n++)
	{
		for (k=0;k<(knode-1);k++)
		{
			msho.iconf[(knode-1)*n+k] = (int)out.trifacelist[(knode-1)*n+k];
		}
		msho.bdflagf[n] = (int)out.trifacemarkerlist[n];
	}
	return;
}

void interpsol(double *sn,int nnode)
{
	int n,el,i,j,inck;
	double v,vk;
	int knode = msh.knode;
	int dim = msh.dim;
	//int nnode = out.numberofpoints;
	double *newu = (double *)malloc(sizeof(double)*nnode);
	double *beta = (double *)malloc(sizeof(double)*knode);
	double **ik = (double **)malloc(sizeof(double *)*knode);
	int count = 0; // using for debugging

	printf("TETGEN: Interpolating to new mesh\n");

	for (n=0;n<nnode;n++) newu[n] = 1.0;

	for (n=0;n<nnode;n++)
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
				ik[i] = &(sn[n*dim]); // only pointer for x coord
				//printf("pk[%d] = [%.5f %.5f %.5f];\n",i,*ik[i],*(ik[i]+1),*(ik[i]+2));
				vk = volume(ik[0],ik[1],ik[2],ik[3]);
				beta[i] = vk/v;
				// dodgy rounding - round to 10 d.p since the vtk is to 6 d.p
				// there was an issue with indentifying the sign of beta
				//double tmp = (long int)(beta[i]*1e10 + .5);
				//beta[i] = tmp / 1e10;

				// checks sign and magnitude
				if (beta[i] < 0.0 && fabs(beta[i]) > 1.e-10)
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
					newu[n] += fabs(beta[i])*spstiff.sol[msh.icon[knode*el+i]];
				}
				count++;
				//if (count==5) printf("	newu[%d] = %.5f\n",n,newu[n]);
			}
		}
	}
	spstiff.sol = (double*)realloc(spstiff.sol,sizeof(double)*nnode);
	for (n=0;n<nnode;n++)
	{
		spstiff.sol[n] = newu[n];
	}

	free(newu); free(beta); free(ik);
	return;
}

void interplocal(double *sn,int nnode)
// using the localisation algorithm
{
	int el,eni,en,found,i,j,n;
	double v,vk;
	int knode = msh.knode;
	int dim = msh.dim;
	int *neigh = msh.neigh;
	 double *newu = (double *)malloc(sizeof(double)*nnode);
        double *beta = (double *)malloc(sizeof(double)*knode);
	double **ik = (double **)malloc(sizeof(double *)*knode);

	printf("TETGEN: Interpolating to new mesh\n");

	// useful for debugging
	for (n=0;n<nnode;n++) newu[n] = 1.0;

	for (n=0;n<nnode;n++)
	//for (n=0;n<1001;n++)
	//for (n=0;n<1;n++)
	{
		eni = 0; // starting element for now
		found = 0;
		//printf("sn = (%.5f, %.5f, %.5f)\n",sn[n*dim],sn[n*dim+1],sn[n*dim+2]);

		while (found != 1)
		{
			en = eni;
			/*
                        printf("p =\n");
                        printf("[%.5f %.5f %.5f;\n",msh.s[msh.icon[knode*en]*dim],msh.s[msh.icon[knode*en]*dim+1],msh.s[msh.icon[knode*en]*dim+2]);
                        printf("%.5f %.5f %.5f;\n",msh.s[msh.icon[knode*en+1]*dim],msh.s[msh.icon[knode*en+1]*dim+1],msh.s[msh.icon[knode*en+1]*dim+2]);
                        printf("%.5f %.5f %.5f;\n",msh.s[msh.icon[knode*en+2]*dim],msh.s[msh.icon[knode*en+2]*dim+1],msh.s[msh.icon[knode*en+2]*dim+2]);
                        printf("%.5f %.5f %.5f];\n",msh.s[msh.icon[knode*en+3]*dim],msh.s[msh.icon[knode*en+3]*dim+1],msh.s[msh.icon[knode*en+3]*dim+2]);
			*/
			v = evolume(en);
			found = 1;
			for (i=0;i<knode;i++)
			{
				for (j=0;j<knode;j++)
				{
					ik[j] = &(msh.s[msh.icon[knode*en+j]*dim]);
				}
				ik[i] = &(sn[n*dim]);
				//printf("pk[%d] = [%.5f %.5f %.5f];\n",i,*ik[i],*(ik[i]+1),*(ik[i]+2));
				vk = volume(ik[0],ik[1],ik[2],ik[3]);
				beta[i] = vk/v;
				//printf("beta[%d] = %.5f\n",i,beta[i]);
				eni = neigh[knode*en+i];
				//printf("en = %d\n",en+1);
				//printf("eni = %d\n",eni+1);

				// some dodgy round, since vtks are only written to 6 d.p.
				if (beta[i]<0. && fabs(beta[i])>1.e-10 && eni>-1 && en!=eni)
				{
					found = 0;
					break;			
				}
				found = 1;
			}
		}
		if (found == 1)
		{
			newu[n] = 0.;
			// interpolation based upon volume ratio
			for (i=0;i<knode;i++)
			{
				newu[n] += fabs(beta[i])*spstiff.sol[msh.icon[knode*en+i]];
			}
			//printf("newu[%d] = %.5f\n",n,newu[n]);
		}
	}
        spstiff.sol = (double*)realloc(spstiff.sol,sizeof(double)*nnode);
        for (n=0;n<nnode;n++)
        {
                spstiff.sol[n] = newu[n];
        }

	free(newu); free(beta); free(ik);
	return;
}
