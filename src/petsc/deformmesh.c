#include "petsc.h"

static void translateregion(double,double *,int *);
static void rotateregion(double,double *,double *,int *);
static void bulgez(double,double *,int *,double,double);

// radial basis interpolation functions
static void RBFinterp(double *,int *);
static void radialsolve(list_node **,double *,double *,double *,int *,int,int,int);
static void cpvec(double *,double *,double *,int *,int,int);
static double phifunction(double);
void poplist(list_node **,double,int);

void deformmesh(int deftype,double dx)
{
	int nnode = msh.nnode;
	int dim = msh.dim;
	int nel = msh.nel;
	int knode = msh.knode;
	int *icon = msh.icon;
	double *s = msh.s;
	int el,n,i;
	int *reg;
	double *sn;

	reg = (int*)calloc(nnode,sizeof(int)); // nodes in the region of deformation
	sn = (double*)malloc(sizeof(double)*nnode*dim); // new mesh

	// be careful with overlapping nodes
	// finding all interior nodes
	for (el=0;el<nel;el++)
	{
		if (msh.region[el]==1)
		{
			for (i=0;i<knode;i++) reg[icon[el*knode+i]] = 1;
		}
	}

	// populate new mesh	
	for (n=0;n<nnode;n++)
	{
		for (i=0;i<dim;i++)
		{
			sn[n*dim+i] = s[n*dim+i];
		}
	}

	switch (deftype)
	{
		case 0:
			translateregion(dx,sn,reg);
			break;
		case 1:
			double xr[3] = {0.5, 0.5, 0.5};
			double rotangr = dx*(M_PI/180.);
			rotateregion(rotangr,xr,sn,reg);
			break;
		case 2:
			bulgez(dx,sn,reg,0.35,0.65);
			break;
	}
	RBFinterp(sn,reg);

	// write over the old mesh
        for (n=0;n<nnode;n++) 
	{
		for (i=0;i<dim;i++) msh.s[n*dim+i] = sn[n*dim+i];
	}

	free(sn);
	free(reg);
	return;
}


static void translateregion(double dx,double *sn,int *reg)
{
	int nnode = msh.nnode;
	int dim = msh.dim;
	double *s = msh.s;
	int n;
	
	// define deformation here
	// shift interior nodes
	for (n=0;n<nnode;n++)
	{
		if (reg[n]==1) 
		{
			sn[n*dim] = s[n*dim] + dx*(s[n*dim+1]-0.35)/0.3;
			//printf("s[%d] = %.5f\n",n,sn[n*dim]);
		}
	}

	return;
}

static void rotateregion(double theta,double xr[3],double *sn,int *reg)
// rotation angle, center of rotation
{
	int nnode = msh.nnode;
	int dim = msh.dim;
	double *s = msh.s;
	int n;
	double yp,zp;
	
	printf("Rotate mesh by %.2f rad...\n",theta);

	for (n=0;n<nnode;n++)
	{
		if (reg[n]==1)
		{
			yp = s[n*dim+1]-xr[1];
			zp = s[n*dim+2]-xr[2];
			sn[n*dim] = s[n*dim];
			sn[n*dim+1] = (cos(theta)*yp)-(sin(theta)*zp)+xr[1];
			sn[n*dim+2] = (sin(theta)*yp)+(cos(theta)*zp)+xr[2];
		}
	}

	return;
}

static void bulgez(double scale,double *sn,int *reg,double z1,double z2)
// buldge the interior faces in the z direction
{
	int nnode = msh.nnode;
	int dim = msh.dim;
	double *s = msh.s;
	int n;

	for (n=0;n<nnode;n++)
	{
		if (reg[n]==1)
		{
			sn[n*dim+2] += scale*(z1-s[n*dim])*(z2-s[n*dim])*(z1-s[n*dim+1])*(z2-s[n*dim+1]);
		}
	}

	return;
}

static void RBFinterp(double *sn,int *reg)
// radial interpolation of the unknown mesh nodes
{
	int nnode = msh.nnode;
	int dim = msh.dim;
	double *s = msh.s;
	int *bd = msh.bdflag;
	int n,i,j,k;
	
	int nb = 0;
	int *bi = (int*)calloc(nnode,sizeof(int));
	int *b = (int*)calloc(nnode,sizeof(int));
	for (n=0;n<nnode;n++)
	{
		if (reg[n]==1 || bd[n]<0) // known nodes (inner deform and boundary)
		{
			b[n] = nb;
			bi[nb] = n;
			nb++;
		}
		else
		{
			//printf("(%.5e, %.5e, %.5e), n = %d\n",s[n*dim],s[n*dim+1],s[n*dim+2],n);
			b[n] = -1;
		}
	}
	bi = (int *)realloc(bi,sizeof(int)*nb);
	int neq = nb+dim+1;
	//printf("neq = %d\n",neq);

	// using linked list
	// only populate the strict upper diag, then use a symmetric array in petsc
	list_node **alist = (list_node **)calloc(neq*dim,sizeof(list_node *));
	int *nnz = (int *)calloc(neq*dim,sizeof(int)); // non-zero in the upper diag and diag
	list_node *current;
	double phi;
	int nz = 0;
	int d,id;
	// compromise between vector of nb*nb and calling phifunc multiple times
	double *phij = (double *)calloc(nb,sizeof(double));

	/*
	for (d=0;d<dim;d++)
	{
	        for (i=0;i<nb;i++)
	        {
			id = i+d*neq;
			alist[id] = (list_node *)malloc(sizeof(list_node));
			alist[id]->na = id;
			alist[id]->val = 1.; // diagonal is always one
			alist[id]->next = NULL;
			current = alist[id];
			nnz[id] = 1;
			//for (j=i+1;j<nb;j++) // only strict upper triangular mat
	                //{
			//	phi = phifunction(&s[bi[i]*dim],&s[bi[j]*dim],dim);
			//	if (phi != 0.0)
			//	{
			//		poplist(&current,phi,j+d*neq);
			//		nnz[id]++;
			//	}
			//}

			phifuncvec(phij,&s[bi[i]*dim],s,&bi[i+1],(nb-i-1),dim);
			for (j=0;j<(nb-i-1);j++) // only strict upper triangular mat
	                {
				if (phij[j] != 0.0)
				{
					poplist(&current,phij[j],j+i+1+d*neq);
					nnz[id]++;
				}
			}			
			// P block
			poplist(&current,1.,nb+d*neq);
			nnz[id]++;
			for (j=0;j<dim;j++)
			{
				poplist(&current,s[bi[i]*dim+j],nb+j+1+d*neq);
				nnz[id]++;
			}
			nz += nnz[id];
		}
		for (i=nb;i<neq;i++) 
		{
			id = i+d*neq;
			// change the zero block to 1e-10 identity
			alist[id] = (list_node *)malloc(sizeof(list_node));
	                alist[id]->na = id;
	                alist[id]->val = 1.e-10; // scaling the diagonal
	                alist[id]->next = NULL;
			nnz[id] = 1;
			nz += nnz[id];
		}
	}
	*/
	for (i=0;i<nb;i++)
	{
		cpvec(phij,&s[bi[i]*dim],s,&bi[i+1],(nb-i-1),dim);
		for (d=0;d<dim;d++)
	        {
			id = i+d*neq;
			alist[id] = (list_node *)malloc(sizeof(list_node));
			alist[id]->na = id;
			alist[id]->val = 1.; // diagonal is always one
			alist[id]->next = NULL;
			current = alist[id];
			nnz[id] = 1;
			for (j=0;j<(nb-i-1);j++) // only strict upper triangular mat
	                {
				if (phij[j] != 0.0)
				{
					poplist(&current,phij[j],j+i+1+d*neq);
					nnz[id]++;
				}
			}			
			// P block
			poplist(&current,1.,nb+d*neq);
			nnz[id]++;
			for (j=0;j<dim;j++)
			{
				poplist(&current,s[bi[i]*dim+j],nb+j+1+d*neq);
				nnz[id]++;
			}
			nz += nnz[id];
		}
	}
	for (i=nb;i<neq;i++) 
	{
		for (d=0;d<dim;d++)
		{
			id = i+d*neq;
			// change the zero block to 1e-10 identity
			alist[id] = (list_node *)malloc(sizeof(list_node));
	                alist[id]->na = id;
	                alist[id]->val = 1.e-10; // scaling the diagonal
	                alist[id]->next = NULL;
			nnz[id] = 1;
			nz += nnz[id];
		}
	}

	//writelist(alist,nz,neq*dim,"alist.dat");
	
	double *f = (double *)calloc(neq*dim,sizeof(double));
	double *alpha = (double *)calloc(nb*dim,sizeof(double));
	double *beta = (double *)calloc((dim+1)*dim,sizeof(double));
	for (d=0;d<dim;d++)
	{
		for (j=0;j<nb;j++) f[j+d*neq] = sn[bi[j]*dim+d];
		for (j=nb;j<neq;j++) f[j+d*neq] = 0.;
		//for (j=0;j<neq;j++) printf("f[%d] = %.5e\n",j+d*neq,f[j+d*neq]);
	}
	radialsolve(alist,f,alpha,beta,nnz,neq,nb,dim);
	//for (i=0;i<nb;i++) printf("alpha[%d] = %.5e\n",i,alpa[i]);
	//for (i=0;i<(dim+1);i++) printf("beta[%d] = %.5e\n",i,beta[i]);

	double *alp,*bet;
	double tmp,p;
	/*
	for (d=0;d<dim;d++)
	{
		alp = &alpha[d*nb];
		bet = &beta[d*(dim+1)];
		for (n=0;n<nnode;n++)
		{
			if (b[n]==-1)
			{
				//tmp = 0.;
				//for (k=0;k<nb;k++)
				//{
				//	tmp += alp[k]*phifunction(&s[n*dim],&s[bi[k]*dim],dim);
				//}
				phifuncvec(phij,&s[n*dim],s,bi,nb,dim);
				tmp = 0.; 
				for (k=0;k<nb;k++) tmp += alp[k]*phij[k];
				p = bet[0] + (s[n*dim]*bet[1]) + (s[n*dim+1]*bet[2]) + (s[n*dim+2]*bet[3]);
				sn[n*dim+d] = tmp + p;
				if (sn[n*dim+d]<0. || sn[n*dim+d]>1.) printf("Error: node leakage, new location outside of the domain.\n");
			}
			//printf("s[%d] = %.5f\n",n,sn[n*dim+i]);
		}
	}
	*/
	for (n=0;n<nnode;n++)
	{
		if (b[n]==-1)
		{
			cpvec(phij,&s[n*dim],s,bi,nb,dim);
			for (d=0;d<dim;d++)
			{
				alp = &alpha[d*nb];
				bet = &beta[d*(dim+1)];
				tmp = 0.; 
				for (k=0;k<nb;k++) tmp += alp[k]*phij[k];
				p = bet[0] + (s[n*dim]*bet[1]) + (s[n*dim+1]*bet[2]) + (s[n*dim+2]*bet[3]);
				sn[n*dim+d] = tmp + p;
                                if (sn[n*dim+d]<0. || sn[n*dim+d]>1.) printf("Error: node leakage, new location outside of the domain.\n");
			}
		}
	}

	free(phij);
	free(f); free(alpha); free(beta);
	free(bi); free(b); free(nnz);
	for (i=0;i<neq*dim;i++)
	{
		listfree(alist[i]);
		alist[i] = NULL;
	}
	free(alist);

	// interpolate solution?

	return;
}

static void cpvec(double *phi,double *xi,double *xs,int *bi,int n,int dim)
// rbf compact support vectorised
{
	int i,j;
	// the larger the value, the more accurate and more dense the solution
	//double r = 0.2; // works well
	double r = 2.5*0.3; // 2.5 times characteristic length
	double x;

	for (i=0;i<n;i++)
	{
		x = 0.;
		for (j=0;j<dim;j++) x += (*(xi+j)-xs[*(bi+i)*dim+j]) * (*(xi+j)-xs[*(bi+i)*dim+j]);
		x = sqrt(x);
		x /= r;
		phi[i] = 0.;
		if (x <= 1.) phi[i] = phifunction(x);
	}
	return;
}

static double phifunction(double x)
// compact support basis functions
{
	double phi;
	//phi = (1.-eps)*(1.-eps); // CP c0
	//phi = pow(1.-eps,4.)*(4.*eps+1.); // c2
	//phi = pow(1.-eps,6.)*(((35./3.)*pow(eps,2.))+6.*eps+1.); // c6
	phi = 1. + (80./3.)*(x*x) - 40.*(x*x*x) + 15.*(x*x*x*x) - (8./3.)*(x*x*x*x*x) + 20.*(x*x)*log(x);
	return phi;
}

static void radialsolve(list_node **list,double *f,double *alpha,double *beta,int *nnz,int neq,int nb,int dim)
// solves for the radial basis function coefficient alpha and beta
{
	PetscErrorCode ierr;
	Mat A;
	Vec b,x;
	KSP ksp;
	PC pc;
	int *icol = (int *)calloc(neq*dim,sizeof(int));
	double *acol = (double *)calloc(neq*dim,sizeof(double));
	int i,d,ncol;
        list_node *current;

	ierr = MatCreate(PETSC_COMM_WORLD,&A);
        ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,neq*dim,neq*dim);
	ierr = MatSetType(A,MATSEQSBAIJ);
	//ierr = MatSetType(A,"sbaij");
	ierr = MatSeqSBAIJSetPreallocation(A,1,PETSC_DEFAULT,nnz); // use array of nnz row allocation

	ierr = VecCreate(PETSC_COMM_WORLD,&x);
        ierr = VecSetSizes(x,PETSC_DECIDE,neq*dim);
        ierr = VecSetFromOptions(x);
        ierr = VecDuplicate(x,&b);

	for (i=0;i<(neq*dim);i++)
	{
		current = list[i];
		ncol = 0;
		while (current!=NULL)
		{
			icol[ncol] = current->na;
			acol[ncol] = current->val;
			ncol++;
			current = current->next;
		}
		ierr = MatSetValues(A,1,&i,ncol,icol,acol,INSERT_VALUES);
		ierr = VecSetValues(b,1,&i,&f[i],INSERT_VALUES);
	}

        ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
        ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
        ierr = VecAssemblyBegin(b);
        ierr = VecAssemblyEnd(b);
	
	// simple cholesky direct solve
	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);
        ierr = KSPSetOperators(ksp,A,A);
	ierr = KSPGetPC(ksp,&pc);
        //ierr = KSPSetType(ksp,KSPPREONLY);
        //ierr = PCSetType(pc,PCCHOLESKY);
        ierr = KSPSetType(ksp,KSPMINRES);
	ierr = PCSetType(pc,PCICC);
	ierr = KSPSetTolerances(ksp,1.e-12,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);

        ierr = KSPSetUp(ksp);
        ierr = KSPSolve(ksp,b,x);

        PetscScalar *get;
        ierr = VecGetArray(x,&get);
	for (d=0;d<dim;d++)
	{
	        for (i=0;i<nb;i++)
		{
			alpha[i+d*nb] = get[i+d*neq];
		}
	        for (i=0;i<(dim+1);i++) beta[i+d*(dim+1)] = get[nb+i+d*neq];
	}


	free(icol); free(acol);
	KSPDestroy(&ksp); MatDestroy(&A); VecDestroy(&b); VecDestroy(&x);
	//PCDestroy(&pc); // no need since part of ksp
	return;
}

