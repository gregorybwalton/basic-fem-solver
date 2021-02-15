#include "petsc.h"

void pfin(PETSC_STRUC *pets)
{
        PetscErrorCode ierr;

	// Clear workspace
	ierr = VecDestroy(&pets->x);
	ierr = VecDestroy(&pets->b);
	ierr = MatDestroy(&pets->A);
	ierr = KSPDestroy(&pets->ksp);
	
        ierr = PetscFinalize();

        printf("Petsc finialised\n");
	
	return;
}
