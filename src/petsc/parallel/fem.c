static char help[] = "Simple 3D FEM solver.\n\n";

#include "../../fem.h"

// Define functions here
void readmesh(char *,char *);
void spmatrices(void);
void ptsolve(int,char **,char *,bool);

int main (int argc, char *argv[])
{
	soltype = 1; // just some boolean which change the problem type
	readmesh("mesh_6","/usr/not-backed-up/phd/fem-codes/fem-working/data/3d/");
	
	spmatrices(); // generates the matrices in spare form
	ptsolve(argc,argv,help,true); // Solve sparse matrix using PETSc, working in serial


	return 0;
}
