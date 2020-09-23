#include "../fem.h"

// Define functions here
void readmesh(char *,char *);
void vtkoutput(double *);
void resoutput(double *);
double calcerror(double *,int);
void spmatrices(void);
void spsolve(struct sysmat);
void printspmatrices(struct sysmat, double *);
void jacobi(struct sysmat,double *);
void ichol(struct sysmat,double *);
char* meshselect(int,char **);

int main (int nargin, char *argsin[])
{
	soltype = 1; // just some boolean which change the problem type
	readmesh(meshselect(nargin,argsin),"../../data/3d/");
	
	spmatrices(); // generates the matrices in spare form
	//printspmatrices(spstiff,load);
	//jacobi(spstiff,load);
	//ichol(spstiff,load);
	spsolve(spstiff); // A sparse matrix version of solve

	vtkoutput(spstiff.sol);
	resoutput(spstiff.sol);
	//calcerror(spstiff.sol,0);

	return 0;
}

char* meshselect(argc,argv)
// checks input for "mesh_*"
int argc;
char** argv;
{
        int i;
        char tstr[10];

        for (i=0;i<argc;i++)
        {
                strncpy(tstr,argv[i],4);
                if (strcmp(tstr,"mesh")==0)
                {
                        return argv[i];
                }
        }
        //return "mesh_2";
        return "mesh_0";
}
