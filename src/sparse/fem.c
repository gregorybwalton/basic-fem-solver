#include "../fem.h"

// Define functions here
void readmesh(char *,char *);
void vtkoutput(double *);
void resoutput(double *);
void calcerror(double *);
void spmatrices(void);
void spsolve(struct sysmat);
void printspmatrices(struct sysmat, double *);
void jacobi(struct sysmat,double *);
void ichol(struct sysmat,double *);


int main ()
{
	soltype = 1; // just some boolean which change the problem type
	//readmesh("inter.1","../../data/region/");
	readmesh("mesh_8","../../data/3d/");
	
	spmatrices(); // generates the matrices in spare form
	//printspmatrices(spstiff,load);
	//jacobi(spstiff,load);
	//ichol(spstiff,load);
	spsolve(spstiff); // A sparse matrix version of solve

	vtkoutput(spstiff.sol);
	resoutput(spstiff.sol);
	calcerror(spstiff.sol);

	return 0;
}
