#include "petsc.h"

void translateregion(double,double *,int *);
void rotateregion(double,double *,double *,int *);
void bulgez(double,double *,int *,double,double);

// radial basis interpolation functions
void radialinterp(double *,int *);
void radialsolve(list_node **,double *,double *,double *,int *,int,int,int);
double phifunction(double *,double *,int);
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
	radialinterp(sn,reg);

	// write over the old mesh
        for (n=0;n<nnode;n++) 
	{
		for (i=0;i<dim;i++) msh.s[n*dim+i] = sn[n*dim+i];
	}

	free(sn);
	free(reg);
	return;
}


void translateregion(double dx,double *sn,int *reg)
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

void rotateregion(double theta,double xr[3],double *sn,int *reg)
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

void bulgez(double scale,double *sn,int *reg,double z1,double z2)
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

void radialinterp(double *sn,int *reg)
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
	list_node **alist = (list_node **)calloc(neq,sizeof(list_node *));
	list_node *current;
	double phi;
	int *nnz = (int *)calloc(neq,sizeof(int)); // non-zero in the upper diag and diag
	int nz = 0;
        for (i=0;i<nb;i++)
        {
		alist[i] = (list_node *)malloc(sizeof(list_node));
		alist[i]->na = i;
		alist[i]->val = 1.; // diagonal is always one
		alist[i]->next = NULL;
		current = alist[i];
		nnz[i] = 1;
		for (j=i+1;j<nb;j++) // only strict upper triangular mat
                {
			phi = phifunction(&s[bi[i]*dim],&s[bi[j]*dim],dim);
			if (phi != 0.0)
			{
				poplist(&current,phi,j);
				nnz[i]++;
			}
		}
		// P block
		poplist(&current,1.,nb);
		nnz[i]++;
		for (j=0;j<dim;j++)
		{
			poplist(&current,s[bi[i]*dim+j],nb+j+1);
			nnz[i]++;
		}
		nz = nz + nnz[i];
	}
	for (i=nb;i<neq;i++) 
	{
		// change the zero block to 1e-10 identity
		alist[i] = (list_node *)malloc(sizeof(list_node));
                alist[i]->na = i;
                alist[i]->val = 1.e-10; // diagonal is always one
                alist[i]->next = NULL;
		nnz[i] = 1;
		//nnz[i]++;
		//nnz[i] = 0; // last 4 rows are empty
		nz = nz + nnz[i];
	}

	//writelist(alist,nz,neq,"alist.dat");

	double *f = (double *)calloc(neq,sizeof(double));
	double *alpha = (double *)calloc(nb,sizeof(double));
	double *beta = (double *)calloc(dim+1,sizeof(double));
	double tmp;
	double p;
	for (i=0;i<dim;i++)
	{
		// displacements at the boundaries
		for (j=0;j<nb;j++) f[j] = sn[bi[j]*dim+i];
		for (j=nb;j<neq;j++) f[j] = 0.;

		// solve the system to find alphas and betas
		radialsolve(alist,f,alpha,beta,nnz,neq,nb,dim);
		//for (j=0;j<nb;j++) printf("alpha[%d] = %.5e\n",j,alpha[j]);
		//for (j=0;j<dim+1;j++) printf("beta[%d] = %.5e\n",j,beta[j]);

		for (n=0;n<nnode;n++)
		{
			if (b[n]==-1)
			{
				tmp = 0.;
				for (k=0;k<nb;k++)
				{
					//tmp += alpha[k]*phifunction(&sn[n*dim],&sn[bi[k]*dim],dim);
					tmp += alpha[k]*phifunction(&s[n*dim],&s[bi[k]*dim],dim);
				}
				p = beta[0] + (s[n*dim]*beta[1]) + (s[n*dim+1]*beta[2]) + (s[n*dim+2]*beta[3]);
				//sn[n*dim+i] = tmp + beta[0] + (sn[n*dim]*beta[1]) + (sn[n*dim+1]*beta[2]) + (sn[n*dim+2]*beta[3]);
				sn[n*dim+i] = tmp + p;
				if (n==1162 && i==1)
				{
					printf("tmp[%d] = %.6f\n",n,tmp);
					printf("p[%d] = %.6f\n",n,p);
				}
				if (sn[n*dim+i]<0. || sn[n*dim+i]>1.) printf("Error: node leakage, new location outside of the domain.\n");
			}
			//printf("s[%d] = %.5f\n",n,sn[n*dim+i]);
		}
	}
	free(f); free(alpha); free(beta);
	free(bi); free(b); free(nnz);
	for (i=0;i<neq;i++)
	{
		listfree(alist[i]);
		alist[i] = NULL;
	}
	free(alist);
	// interpolate solution?

	return;
}

double phifunction(double *xi,double *xj,int dim)
// definition of the radial basis function
{
	double phi;
	int i;
	// rule of thumb: smallest as possible
	double r = 0.05; // radius of influence
	double x = 0.;
	for (i=0;i<dim;i++) x += (*(xi+i)-*(xj+i))*(*(xi+i)-*(xj+i));
	x = sqrt(x);
	x /= r;
	if (x > 1){
		phi = 0.;
	} else {
		phi = (1.-x)*(1.-x); // radial basis function compact support
		//phi = pow(1.-x,4.)*(4.*x+1.); // c2
		//phi = pow(1.-x,6.)*(((35./3.)*pow(x,2.))+6.*x+1.); // c6
	}

	return phi;
}

void poplist(list_node **curr,double cval,int cna)
// used to build and populate a sparse linked list
// worth moving
{
        (*curr)->next = (list_node *)malloc(sizeof(list_node));
	(*curr)->next->na = cna;
	(*curr)->next->val = cval;
	(*curr)->next->next = NULL;
	(*curr) = (*curr)->next;
	return;
}



void radialsolve(list_node **list,double *f,double *alpha,double *beta,int *nnz,int neq,int nb,int dim)
// solves for the radial basis function coefficient alpha and beta
{
	PetscErrorCode ierr;
	Mat A;
	Vec b,x;
	KSP ksp;
	PC pc;
	int *icol = (int *)calloc(neq,sizeof(int));
	double *acol = (double *)calloc(neq,sizeof(double));
	int i,ncol;
        list_node *current;

	ierr = MatCreate(PETSC_COMM_WORLD,&A);
        ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,neq,neq);
	ierr = MatSetType(A,MATSEQSBAIJ);
	//ierr = MatSetType(A,"sbaij");
	ierr = MatSeqSBAIJSetPreallocation(A,1,PETSC_DEFAULT,nnz); // use array of nnz row allocation

	ierr = VecCreate(PETSC_COMM_WORLD,&x);
        ierr = VecSetSizes(x,PETSC_DECIDE,neq);
        ierr = VecSetFromOptions(x);
        ierr = VecDuplicate(x,&b);

	for (i=0;i<neq;i++)
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
	ierr = KSPSetTolerances(ksp,1.e-2,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);

        ierr = KSPSetUp(ksp);
        ierr = KSPSolve(ksp,b,x);

        PetscScalar *get;
        ierr = VecGetArray(x,&get);
        for (i=0;i<nb;i++) alpha[i] = get[i];
        for (i=0;i<(dim+1);i++) beta[i] = get[nb+i];


	free(icol); free(acol);
	KSPDestroy(&ksp); MatDestroy(&A); VecDestroy(&b); VecDestroy(&x);
	//PCDestroy(&pc); // no need since part of ksp
	return;
}

