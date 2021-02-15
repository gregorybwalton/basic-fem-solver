#include "petsc.h"

void pinit(PETSC_STRUC *pets,int *argc,char ***argv,PetscMPIInt *size)
{
        PetscErrorCode ierr;
        PetscMPIInt rank;
	int *nnz = spstiff.nnzeros;
	PetscInt neq = msh.neq; //number of unknowns
	
        ierr = PetscInitialize(argc,argv,(char*) 0,(char*) 0);
        ierr = MPI_Comm_size(PETSC_COMM_WORLD,size);
        ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

        PetscPrintf(PETSC_COMM_WORLD,"Petsc initialised, with %d processes\n",*size);
        ierr = PetscMemorySetGetMaximumUsage(); // tells petsc to monitor memory
	
	PetscPrintf(PETSC_COMM_WORLD,"Vec setup...\n");
	ierr = VecCreate(PETSC_COMM_WORLD,&pets->x);
	ierr = VecSetSizes(pets->x,PETSC_DECIDE,neq);
	ierr = VecSetType(pets->x,"standard");
	ierr = VecSetFromOptions(pets->x);
	ierr = VecDuplicate(pets->x,&pets->b);
	
	PetscPrintf(PETSC_COMM_WORLD,"Mat setup...\n");
	// series allocation
	//ierr = MatCreateSeqAIJ(comm,n,n,0,nnz,&A); // compressed row sparse format
	// parallel allocation - See page 95 of the manual
	PetscInt dnz = (neq/(*size));
	PetscInt onz = neq-(neq/(*size));
	ierr = MatCreate(PETSC_COMM_WORLD,&pets->A);
	ierr = MatSetSizes(pets->A,PETSC_DECIDE,PETSC_DECIDE,neq,neq);
	ierr = MatSetFromOptions(pets->A);
	ierr = MatMPIAIJSetPreallocation(pets->A,dnz,NULL,onz,NULL);
	ierr = MatSeqAIJSetPreallocation(pets->A,0,nnz);
	ierr = MatSetUp(pets->A);

	return;
}
