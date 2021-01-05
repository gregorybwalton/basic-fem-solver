static char helpstr[] = "Simple 3D FEM solver.\n\n";

#include "../fem.h"

// Define functions here
void readmesh(char *,char *);
void readvtk(char *,char *);
void vtkoutput(double *);
void resoutput(double *);
void d3n4(void);
void ptsolve(int *,char ***,char *,bool);
char* meshselect(int,char **);
void meshvolume();

void runtetgen();

int main (int nargin, char *argsin[])
{

	soltype = 1; // just some boolean which change the problem type
	//readmesh(meshselect(nargin,argsin),"../../data/3d/");
	readmesh(meshselect(nargin,argsin),"../../data/3d/region/innersmall/");
	//readvtk("test3d","../../data/grummp/");
	//readvtk("inner","../../data/grummp/");

	d3n4(); // generates the matrices in spare form
	ptsolve(&nargin,&argsin,helpstr,false); // Solve sparse matrix using PETSc, working in serial

	//vtkoutput(spstiff.sol);
	//resoutput(spstiff.sol);

	// Adaptivity bit
	meshvolume(); // Calculting the volume requirement
	//meshlengthscale(); // GRUMMP refinement
	runtetgen(); // Running TetGen from code
	vtkoutput(spstiff.sol);

	//free(msh.icon); free(msh.s); free(msh.bdflag); free(msh.region); free(msh.lscon);
	return 0;
}

// checks input for "mesh_*"
char* meshselect(int argc, char** argv)
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
	//return "mesh_9";
	return "mesh_0";
}
