#include "../fem.h"
#include "petscksp.h" //  includes petscmat.h
#include "petscviewer.h"
#include "unistd.h"

PetscErrorCode writebin(Mat*,Vec*,Vec*);
void writeres(char *,int,double,int,double,double);
double calcerror(double *,int);
double droptolselect(int,char **);

void ptsolve(argc,argv,help,wbin)
int *argc;
char ***argv;
char *help;
bool wbin;
{
	//int nz = spstiff.nzeros;
	int *nnz = spstiff.nnzeros;
	int *bdflag = msh.bdflag;
	int nnode = msh.nnode;
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
	double timet,nerr,mem;

	//MPI_Comm comm = MPI_COMM_WORLD;
	Mat A; //linear system matrix
	Vec b,x,u; //petsc rhs vector, solution vector
	KSP ksp; //linear solver
	PC pc; // preconditioner
	PetscReal norm;
	PetscInt its;
	PetscErrorCode ierr;
	PetscMPIInt size, rank;
	PetscScalar neg_one = -1.0;
	PetscScalar one = 1.0;

	ierr = PetscInitialize(argc,argv,(char*) 0,help);
	//ierr = PetscInitialize(NULL,NULL,NULL,help);
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);
	
	PetscPrintf(PETSC_COMM_WORLD,"Vec setup...\n");
	ierr = VecCreate(PETSC_COMM_WORLD,&x);
	ierr = VecSetSizes(x,PETSC_DECIDE,neq);
	ierr = VecSetFromOptions(x);
	ierr = VecDuplicate(x,&b);
	ierr = VecDuplicate(x,&u);
	
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
				//printf("(%d,%d) = %.15f\n",nbd,icol[ncol],vcol[ncol]);
				ncol++;
				current = current->next;
			}
			ierr = MatSetValues(A,1,&nbd,ncol,icol,vcol,INSERT_VALUES);
			ierr = VecSetValues(b,1,&nbd,&spstiff.load[nbd],INSERT_VALUES);
			ierr = VecSetValues(x,1,&nbd,&spstiff.sol[node],INSERT_VALUES);
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
	//ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD);
	//ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);
	
	//check is symmetric
	PetscReal symtol = 0.0;
	PetscBool symflg;
	ierr = MatIsSymmetric(A,symtol,&symflg);
	PetscPrintf(PETSC_COMM_WORLD,"Symmetric = %d\n",symflg);
	
	ierr = VecSet(u,one);
	ierr = MatMult(A,u,b);
	
	if (wbin)
	{
		ierr = writebin(&A,&x,&b);
		PetscPrintf(PETSC_COMM_WORLD,"... printing bin files\n");
		PetscFinalize();
		return;
	}
		
	int soltyp = 3;
	int pctyp = 3;
	// Setup the ksp
	PetscPrintf(PETSC_COMM_WORLD,"Setuping up solver...\n");
	ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);
	ierr = KSPSetOperators(ksp,A,A); // Use the matrix as the preconditioning matrix
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
	ierr = PCFactorSetDropTolerance(pc,droptolselect(*argc,*argv),0.05,dnz);
	ierr = KSPSetTolerances(ksp,1.0e-5,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);

	/*
	PetscReal abstol;
	ierr = KSPGetTolerances(ksp,NULL,&abstol,NULL,NULL);
	printf("%g\n",abstol);
	*/
	
	// set non-zero initial guess here
	// -- Currently zero
	
	// Solve linear system here
	PetscPrintf(PETSC_COMM_WORLD,"Running solver...\n");
	getrusage(RUSAGE_SELF, &usage);
	start = usage.ru_utime;

	ierr = KSPSolve(ksp,b,x);
	ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD); // view solver info

	getrusage(RUSAGE_SELF, &usage);
	end = usage.ru_utime;
	timet = (double)(end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) * 1e-6;
	PetscPrintf(PETSC_COMM_WORLD,"Time taken to solve (CPU time): %.5f s\n", timet);

	mem = (usage.ru_maxrss*0.001); // in MB
	PetscPrintf(PETSC_COMM_WORLD,"Peak RAM Usage: %.2f MB\n", mem);

	// Error check
	ierr = VecAXPY(x,neg_one,u);
	ierr = VecNorm(x,NORM_2,&norm);
	ierr = KSPGetIterationNumber(ksp,&its);
	ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %.15f, Iterations %d\n",(double)norm,its);

	// get solution
	PetscScalar *get;
	ierr = VecGetArray(x,&get);
	for (node=0;node<nnode;node++)
	{
		nbd = bdflag[node];
		if (nbd>=0)
		{
			spstiff.sol[node] = get[nbd];
		}
	}
	ierr = VecRestoreArray(x,&get);

	// Clear workspace
	ierr = VecDestroy(&x);
	ierr = VecDestroy(&b);
	ierr = VecDestroy(&u);
	ierr = MatDestroy(&A);
	ierr = KSPDestroy(&ksp);
	
	PetscPrintf(PETSC_COMM_WORLD,"--------------------\n");
	PetscFinalize();
	
	nerr = calcerror(spstiff.sol,1);
	writeres("./results/study",neq,nerr,its,timet,mem);

	return;
}

PetscErrorCode writebin(A,x,b)
Mat *A;
Vec *x,*b;
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

void writeres(wfil,ukn,norm,its,timet,mem)
char wfil[50];
int ukn,its;
double norm,timet,mem;
{
	FILE *fp;
	if (access(wfil,W_OK) != -1)
	{
		fp = fopen(wfil,"a");
		fprintf(fp,"%i,%i,%g,%g,%g\n",ukn,its,norm,timet,mem);
		fclose(fp);
	}
}

double droptolselect(argc,argv)
// checks input for "-droptol"
int argc;
char** argv;
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
