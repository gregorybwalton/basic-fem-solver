//#include "../fem.h"
#include "petsc.h"

static char helpstr[] = "Simple 3D FEM solver.\n\n";

// Define functions here
void readmesh(char *,char *);
void readvtk(char *,char *);
void vtkoutput(double *,int);
void resoutput(double *);
void d3n4(void);
char* meshselect(int,char **);
void runtetgen();
void deformmesh(int,double);
void writetetgen(char *);
void meshqual(void);

void pinit(int *,char ***,PetscMPIInt *);
void ptsolve(bool,PetscMPIInt);
void pfin();

int main (int argc, char **args)
{
	int i;
	int imax = 16;
	PetscMPIInt psize;
	soltype = 1; // just some boolean which change the problem type
	//readmesh(meshselect(nargin,argsin),"../../data/3d/");
	readmesh(meshselect(argc,args),"../../data/3d/region/innersmall/");

	d3n4(); // generates the matrices in spare form
	pinit(&argc,&args,&psize);
	ptsolve(false,psize); // Solve sparse matrix using PETSc, working in serial

	vtkoutput(spstiff.sol,0);
	meshqual();
	//resoutput(spstiff.sol);

	runtetgen();
	for (i=1;i<=imax;i++)
	{
		printf("\nSTEP %d\n",i);
		//deformmesh(1,.5); // rotation
		deformmesh(2,2.); // bugle deformation
		meshqual();
		//writetetgen("./testoutput/testoutput");
		if ((i%4)==0) {
			runtetgen(); // Running TetGen from code
			meshqual();
		}
		vtkoutput(spstiff.sol,i);
	}
	pfin();

	//meshvolumescale(); // Calculting the volume requirement about the distance to the interior surface
        //meshlengthscale(); // GRUMMP refinement
	//runtetgen(); // Running TetGen from code
	//vtkoutput(spstiff.sol,1);

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
