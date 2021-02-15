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
void writevol();
void meshquality(double *,double *);
void checkqualtetgen(void);

void pinit(PETSC_STRUC *,int *,char ***,PetscMPIInt *);
void ptsolve(PETSC_STRUC *,bool);
void pfin(PETSC_STRUC *);

PETSC_STRUC petsol;

int main (int argc, char **args)
{
	int i;
	int imax = 16;
	//int imax = 90;
	PetscMPIInt psize;
	soltype = 1; // just some boolean which change the problem type
	//readmesh(meshselect(nargin,argsin),"../../data/3d/");
	readmesh(meshselect(argc,args),"../../data/3d/region/innersmall/");
	//readmesh("mesh_4","../../data/3d/region/innersmall/");

	d3n4(); // generates the matrices in spare form
	pinit(&petsol,&argc,&args,&psize);
	ptsolve(&petsol,false); // Solve sparse matrix using PETSc, working in serial

	vtkoutput(spstiff.sol,0);
	//resoutput(spstiff.sol);

	runtetgen();
	for (i=1;i<=imax;i++)
	{
		printf("\nSTEP %d\n",i);
		//deformmesh(0,0.1); // shift deform
		//deformmesh(1,1.); // rotation
		deformmesh(2,5.); // bugle deformation
		
		meshquality(NULL,NULL);
		if ((i%4)==0) {
			runtetgen(); // Running TetGen from code
			meshquality(NULL,NULL);
		}
		
		//checkqualtetgen();
		writevol();
		vtkoutput(spstiff.sol,i);
	}
	pfin(&petsol);

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
	char tstr[4];

	for (i=0;i<argc;i++)
	{
		strncpy(tstr,argv[i],4);
		if (strcmp(tstr,"mesh")==1)
		{
			return argv[i];
		}
	}
	//return "mesh_9";
	return "mesh_0";
}
