#include "fem.h"

double easpectratio(int);

void vtkoutput(double* res,int idx)
{
        int nnode = msh.nnode;
        int nel = msh.nel;
        int knode = msh.knode;
	int dim = msh.dim;
	int *bdflag = msh.bdflag;
        int *icon = msh.icon;
        double *s = msh.s;
	int *region = msh.region;
        int node,el,i;
	double vol,ar;

        char vfnm[50];
	sprintf(vfnm,"%s.vtk.%i",savename(),idx);
	
	printf("Saving: %s\n",vfnm);

        FILE *vtkfil;
        vtkfil = fopen(vfnm,"w");

        int size = nel*(knode+1);
        int cell_type = 10; // 10 is cell type VTK_TETRA, see vtk docs for full detail

        fprintf(vtkfil,"# vtk DataFile Version 2.0\n");
        fprintf(vtkfil,"VTK example\n");
        fprintf(vtkfil,"ASCII\n");
        fprintf(vtkfil,"DATASET UNSTRUCTURED_GRID\n");

        fprintf(vtkfil,"POINTS %d %s\n",nnode,"float");
        for (node=0;node<nnode;node++){
                fprintf(vtkfil,"%2.6f %2.6f %2.6f\n",s[node*dim+0],s[node*dim+1],s[node*dim+2]);
        }
        fprintf(vtkfil,"\n");

        fprintf(vtkfil,"CELLS %d  %d\n",nel,size);
        for (el=0;el<nel;el++){
                fprintf(vtkfil,"%d",knode);
                for (i=0;i<knode;i++){
                        fprintf(vtkfil,"%10d",icon[el*knode+i]);
                }
                fprintf(vtkfil,"\n");
        }
        fprintf(vtkfil,"\n");

	// cell type
        fprintf(vtkfil,"CELL_TYPES %d\n",nel);
        for (el=0;el<nel;el++){
                fprintf(vtkfil,"%d\n",cell_type);
        }
        fprintf(vtkfil,"\n");

	// result data
        fprintf(vtkfil,"POINT_DATA %d\n",nnode);
        fprintf(vtkfil,"SCALARS %s %s %d\n","resultData","float",1);
        fprintf(vtkfil,"LOOKUP_TABLE %s\n","testData");
        for (i=0;i<nnode;i++){
                fprintf(vtkfil,"%2.6f\n",res[i]);
        }
        //fprintf(vtkfil,"\n");

	// boundary flags
        fprintf(vtkfil,"SCALARS %s %s %d\n","boundaryNodes","int",1);
        fprintf(vtkfil,"LOOKUP_TABLE %s\n","default");
        for (i=0;i<nnode;i++){
		if (bdflag[i]==-1) {
	                fprintf(vtkfil,"%2d\n",-1);
		}
		else {
			fprintf(vtkfil,"%2d\n",1);
		}
        }
        fprintf(vtkfil,"\n");


	// print index for each part
	fprintf(vtkfil,"CELL_DATA %d\n",nel);
        fprintf(vtkfil,"SCALARS %s %s %i\n","regions","int",1);
        fprintf(vtkfil,"LOOKUP_TABLE %s\n","mat");
        for (el=0;el<nel;el++)
        {
                fprintf(vtkfil,"%10d\n",region[el]);
	}
        //fprintf(vtkfil,"\n");
        
	// volume of each cell
	fprintf(vtkfil,"SCALARS %s %s %i\n","cellVol","float",1);
	fprintf(vtkfil,"LOOKUP_TABLE %s\n","default");
        for (el=0;el<nel;el++)
	{
		vol = evolume(el);
		fprintf(vtkfil,"%2.6f\n",vol);
	}
	//fprintf(vtkfil,"\n");

	// length ratio of each cell
	fprintf(vtkfil,"SCALARS %s %s %i\n","cellAR","float",1);
	fprintf(vtkfil,"LOOKUP_TABLE %s\n","default");
	for (el=0;el<nel;el++)
	{
		ar = easpectratio(el);
		fprintf(vtkfil,"%2.6f\n",ar);
	}


        fclose(vtkfil);
	return;
}

