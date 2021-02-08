#include "petsc.h"

void pinit(int *argc,char ***argv,PetscMPIInt *size)
{
        PetscErrorCode ierr;
        PetscMPIInt rank;
	
        ierr = PetscInitialize(argc,argv,(char*) 0,(char*) 0);
        ierr = MPI_Comm_size(PETSC_COMM_WORLD,size);
        ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

        PetscPrintf(PETSC_COMM_WORLD,"Petsc initialised, with %d processes\n",*size);

	return;
}
