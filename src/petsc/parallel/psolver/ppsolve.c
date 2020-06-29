//#include "/home/cserv1_a/soc_pg/scgwa/not-backed-up/phd/fem-codes/fem-working/src/fem.h"
#include "petscksp.h"

const char help[] = "Separate PETSc solver, reads in binary file.\n For parallel solving.";

PetscErrorCode readbin(Mat*,Vec*,Vec*);

int main (int argc,char *argv[])
{
PetscErrorCode ierr;
PetscMPIInt size, rank;
Mat A;
Vec x,b,u;
KSP ksp; //linear solver
PC pc; // preconditioner
PetscReal norm;
PetscInt its;
PetscScalar neg_one = -1.0;
PetscScalar one = 1.0;
struct timeval start, end;
struct rusage usage;
double timet;

PetscInitialize(&argc,&argv,(char*) 0,help);
ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);
ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);

PetscPrintf(PETSC_COMM_WORLD,"Loading bin files...\n");
ierr = VecCreate(PETSC_COMM_WORLD,&x);
ierr = VecCreate(PETSC_COMM_WORLD,&b);
ierr = MatCreate(PETSC_COMM_WORLD,&A);
ierr = readbin(&A,&x,&b);
ierr = VecDuplicate(x,&u);

// prints the vecs and mats
//ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD);
//ierr = VecView(x,PETSC_VIEWER_STDOUT_WORLD);
//ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD);

ierr = VecSet(u,one);
ierr = MatMult(A,u,b);

PetscPrintf(PETSC_COMM_WORLD,"Setup solver...\n");
// Setup the ksp
ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);
ierr = KSPSetOperators(ksp,A,A); // Use the matrix as the preconditioning matrix
// Set preconditioner and tolerance
ierr = KSPGetPC(ksp,&pc);
ierr = PCSetType(pc,PCJACOBI);
ierr = PCSetType(pc,PCNONE);
ierr = KSPSetTolerances(ksp,1.0e-5,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);

// set non-zero initial guess here
// -- Currently zero

// Solve linear system here
getrusage(RUSAGE_SELF, &usage);
start = usage.ru_utime;
ierr = KSPSolve(ksp,b,x);
ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD); // view solver info
getrusage(RUSAGE_SELF, &usage);
end = usage.ru_utime;
timet = (double)(end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) * 1e-6;
PetscPrintf(PETSC_COMM_WORLD,"Time taken to solve (CPU time): %.5f s\n", timet);

// Error check
ierr = VecAXPY(x,neg_one,u);
ierr = VecNorm(x,NORM_2,&norm);
ierr = KSPGetIterationNumber(ksp,&its);
ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %.15f, Iterations %d\n",(double)norm,its);

// Clear workspace
ierr = VecDestroy(&x);
ierr = VecDestroy(&b);
ierr = VecDestroy(&u);
ierr = MatDestroy(&A);
ierr = KSPDestroy(&ksp);

PetscFinalize();
return 0;
}

PetscErrorCode readbin(A,x,b)
Mat *A;
Vec *x,*b;
{
PetscViewer vin;
PetscErrorCode ierr;

ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"../x.bin",FILE_MODE_READ,&vin);
ierr = VecLoad(*x,vin);

ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"../b.bin",FILE_MODE_READ,&vin);
ierr = VecLoad(*b,vin);

ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,"../A.bin",FILE_MODE_READ,&vin);
ierr = MatLoad(*A,vin);

PetscViewerDestroy(&vin);
return ierr;


}
