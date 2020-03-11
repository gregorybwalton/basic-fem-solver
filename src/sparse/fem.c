#include "../fem.h"

// Define functions here
void readmesh(char *,char *);
void vtkoutput(double *);
void resoutput(double *);
void calcerror(double *);
void spmatrices(void);
void spsolve(struct spmat,double *,double *);
void printspmatrices(struct spmat, double *);
void jacobi(struct spmat,double *);
void ichol(struct spmat,double *);


int main ()
{
	soltype = 1; // just some boolean which change the problem type
	readmesh("inter.1","../../data/region/");
	//readmesh("mesh_0","../data/3d/");
	
	spmatrices(); // generates the matrices in spare form
	//printspmatrices(spstiff,load);
	jacobi(spstiff,load);
	//ichol(spstiff,load);
	spsolve(spstiff,load,sol); // A sparse matrix version of solve

	vtkoutput(sol);
	resoutput(sol);
	calcerror(sol);

	return 0;
}
