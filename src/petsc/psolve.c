#include "petsc.h"

PetscErrorCode writebin(Mat*,Vec*,Vec*);
void writeres(char *,int,double,int,double,double);
double calcerror(double *,int);
double droptolselect(int,char **);

void ptsolve(bool wbin,PetscMPIInt size)
{
	//int nz = spstiff.nzeros;
	int *nnz = spstiff.nnzeros;
	int *bdflag = msh.bdflag;
	int nnode = msh.nnode;
	int dim = msh.dim;
	int node;
	PetscInt neq = msh.neq; //number of unknowns
	list_node **list = spstiff.head;
	list_node *current;
	int nbd;
	int ncol;
	int icol[neq];
	double vcol[neq];
	struct timeval start, end;
	struct rusage usage;
	double timet,nerr;

	//MPI_Comm comm = MPI_COMM_WORLD;
	Mat A; //linear system matrix
	Vec b,x; //petsc rhs vector, solution vector
	KSP ksp; //linear solver
	PC pc; // preconditioner
	PetscReal norm,rnorm;
	PetscInt its;
	PetscErrorCode ierr;
	PetscScalar neg_one = -1.0;
	PetscScalar one = 1.0;
	PetscReal tol = 1.0e-6;
	PetscInt maxit = 10000;
        PetscReal emax,emin;
	PetscLogDouble mem;

	//ierr = PetscInitialize(argc,argv,(char*) 0,help);
	//ierr = PetscInitialize(NULL,NULL,NULL,help);
	//ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);
	//ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

        ierr = PetscMemorySetGetMaximumUsage(); // tells petsc to monitor memory
	
	PetscPrintf(PETSC_COMM_WORLD,"Vec setup...\n");
	ierr = VecCreate(PETSC_COMM_WORLD,&x);
	ierr = VecSetSizes(x,PETSC_DECIDE,neq);
	ierr = VecSetFromOptions(x);
	ierr = VecDuplicate(x,&b);
	
	PetscPrintf(PETSC_COMM_WORLD,"Mat setup...\n");
	// series allocation
	//ierr = MatCreateSeqAIJ(comm,n,n,0,nnz,&A); // compressed row sparse format
	// parallel allocation - See page 95 of the manual
	PetscInt dnz = (neq/size);
	PetscInt onz = neq-(neq/size);
	ierr = MatCreate(PETSC_COMM_WORLD,&A);
	ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,neq,neq);
	ierr = MatSetFromOptions(A);
	ierr = MatMPIAIJSetPreallocation(A,dnz,NULL,onz,NULL);
	ierr = MatSeqAIJSetPreallocation(A,0,nnz);
	ierr = MatSetUp(A);

	// populate dense matrix, and vectors
	PetscPrintf(PETSC_COMM_WORLD,"Populate mats and vecs...\n");
	for (node=0;node<nnode;node++)
	{
		nbd = bdflag[node];
		if (nbd>=0)
		{
			current = list[nbd];
			ncol = 0;
			while (current != NULL)
			{
				icol[ncol] = current->na;
				vcol[ncol] = current->val;
				ncol++;
				current = current->next;
			}
			ierr = MatSetValues(A,1,&nbd,ncol,icol,vcol,INSERT_VALUES);
			ierr = VecSetValues(b,1,&nbd,&spstiff.load[nbd],INSERT_VALUES);
			ierr = VecSetValues(x,1,&nbd,&spstiff.sol[node],INSERT_VALUES);
			//PetscPrintf(PETSC_COMM_WORLD,"load[%i] = %g\n",nbd,spstiff.load[nbd]);
			//PetscPrintf(PETSC_COMM_WORLD,"sol[%i] = %g\n",nbd,spstiff.sol[node]);
		}
	}
	
	ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
	ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
	
	ierr = VecAssemblyBegin(b);
	ierr = VecAssemblyEnd(b);
	ierr = VecAssemblyBegin(x);
	ierr = VecAssemblyEnd(x);
	
	// prints the vecs and mats
	//ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD);
	//ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);
	//ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD);
	
	//check is symmetric
	PetscReal symtol = 0.0;
	PetscBool symflg;
	ierr = MatIsSymmetric(A,symtol,&symflg);
	PetscPrintf(PETSC_COMM_WORLD,"Symmetric = %d\n",symflg);

	if (wbin)
	{
		ierr = writebin(&A,&x,&b);
		PetscPrintf(PETSC_COMM_WORLD,"... printing bin files\n");
		PetscFinalize();
		return;
	}
		
	// Setup the solver
	PetscPrintf(PETSC_COMM_WORLD,"Setuping up solver...\n");
	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);
	ierr = KSPSetOperators(ksp,A,A); // Use the matrix as the preconditioning matrix
	//ierr = PCFactorSetDropTolerance(pc,droptolselect(*argc,*argv),0.05,dnz);
	ierr = KSPSetTolerances(ksp,tol,PETSC_DEFAULT,PETSC_DEFAULT,maxit);
	// set non-zero initial guess here
	// -- Currently zero

	int soltyp = 3;
	int pctyp = 1;
	switch (soltyp)
	{
		case 0:
			ierr = KSPSetType(ksp,KSPPREONLY); // Applied only the preconditioner once
			// May only working with PCLU and PCCHOLESKY
			break;
		case 1:
			ierr = KSPSetType(ksp,KSPCG); // Set solver
			break;
		case 2:
			ierr = KSPSetType(ksp,KSPBCGS);
			break;	
		case 3:
			ierr = KSPSetType(ksp,KSPGMRES);
			break;
		case 4:
			ierr = KSPSetType(ksp,KSPMINRES);
			break;
	}
	

	// Set preconditioner and tolerance
	ierr = KSPGetPC(ksp,&pc);
	switch (pctyp)
	{
		case 0:
			ierr = PCSetType(pc,PCNONE);
			break;
		case 1:
			ierr = PCSetType(pc,PCJACOBI);
			break;
		case 2:
			if (soltyp==0)
			{
				ierr = PCSetType(pc,PCLU); // complete LU
			}
			else
			{
				ierr = PCSetType(pc,PCILU); // incomplete LU
			}
			break;
		case 3:
			if (soltyp==0)
			{
				ierr = PCSetType(pc,PCCHOLESKY);
			}
			else
			{
				ierr = PCSetType(pc,PCICC);
			}
			break;
		case 4:
			ierr = PCSetType(pc,PCGAMG);
			ierr = PCGAMGSetType(pc,PCGAMGAGG);
			break;
	}	

	// for condition number
        ierr = KSPSetComputeSingularValues(ksp,PETSC_TRUE);
        ierr = KSPSetUp(ksp);
	
	// Solve linear system here
	PetscPrintf(PETSC_COMM_WORLD,"Running solver...\n");
	getrusage(RUSAGE_SELF, &usage);
	start = usage.ru_utime;
	ierr = KSPSolve(ksp,b,x);
	getrusage(RUSAGE_SELF, &usage);
	end = usage.ru_utime;
	timet = (double)(end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) * 1e-6;

	ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD); // view solver info

        // Calculate condition number
        ierr = KSPComputeExtremeSingularValues(ksp,&emax,&emin);
        double acond = emax/emin;
        
        ierr = KSPGetIterationNumber(ksp,&its);
        ierr = KSPGetResidualNorm(ksp,&rnorm);
        ierr = PetscMemoryGetMaximumUsage(&mem);

        PetscPrintf(PETSC_COMM_WORLD,"Residual norm = %14.6e\nIterations = %i\n",(double)rnorm,its);
        PetscPrintf(PETSC_COMM_WORLD,"Condition number system = %.2f\n",acond);
        PetscPrintf(PETSC_COMM_WORLD,"Time taken to solve (CPU time, s) = %.5f\n", timet);
        PetscPrintf(PETSC_COMM_WORLD,"Maximum memory usage (MB) = %.2f\n", (double)mem*1.0e-6);

        // Recover solution
	PetscScalar *get;
	ierr = VecGetArray(x,&get);
	for (node=0;node<nnode;node++)
	{
		nbd = bdflag[node];
		if (nbd>=0)
		{
			spstiff.sol[node] = get[nbd];
			//printf("sol(%.5f,%.5f,%.5f) = %.15f\n",msh.s[node*dim+0],msh.s[node*dim+1],msh.s[node*dim+2],spstiff.sol[node]);
		}
	}
	ierr = VecRestoreArray(x,&get);
	
	// Clear workspace
	ierr = VecDestroy(&x);
	ierr = VecDestroy(&b);
	ierr = MatDestroy(&A);
	ierr = KSPDestroy(&ksp);
	
	PetscPrintf(PETSC_COMM_WORLD,"--------------------\n");
	//PetscFinalize();
	
	nerr = calcerror(spstiff.sol,0);
	//writeres("./results/study",neq,nerr,its,timet,mem);

	return;
}

PetscErrorCode writebin(Mat* A,Vec* x,Vec* b)
{
	PetscViewer vout;
	PetscErrorCode ierr;

	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"x.bin",FILE_MODE_WRITE,&vout);
	ierr = VecView(*x,vout);

	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"b.bin",FILE_MODE_WRITE,&vout);
	ierr = VecView(*b,vout);

	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"A.bin",FILE_MODE_WRITE,&vout);
	ierr = MatView(*A,vout);

	PetscViewerDestroy(&vout);
	return ierr;
}

void writeres(char* wfil,int ukn,double norm,int its,double timet,double mem)
{
	FILE *fp;
	if (access(wfil,W_OK) != -1)
	{
		fp = fopen(wfil,"a");
		fprintf(fp,"%i,%i,%g,%g,%g\n",ukn,its,norm,timet,mem);
		fclose(fp);
	}
}

double droptolselect(int argc,char** argv)
// checks input for "-droptol"
{
	int i;
	for (i=0;i<argc;i++)
	{
		if (strcmp(argv[i],"--droptol")==0)
		{
			return strtod(argv[i+1],&argv[i+2]);
		}
	}
	return 1e-3;
}
