#include "../fem.h"
#include "petscksp.h" // includes petscmat.h
#include "petscviewer.h"
#include "petscsys.h"

typedef Mat PETSC_MAT;
typedef Vec PETSC_VEC;

typedef struct
{
        PETSC_MAT A;
        PETSC_VEC x;
        PETSC_VEC b;
        KSP ksp;
        PC pc;
}PETSC_STRUC;
