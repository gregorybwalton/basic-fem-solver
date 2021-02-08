#include "petsc.h"

void pfin()
{
        PetscErrorCode ierr;
        ierr = PetscFinalize();

        printf("Petsc finialised\n");
	
	return;
}
