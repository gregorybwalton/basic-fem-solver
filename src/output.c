#include "fem.h"

char sfnm[50];

char *savename(void)
{
	char tfnm[50];
	char *sptr = sfnm;
	if (soltype==0)
	{
		strcpy(tfnm,"linear");
	}
	else if (soltype==1)
	{
		strcpy(tfnm,"quad");
	}
	//sprintf(sfnm,"%s/%s_%d","./results",tfnm,msh.nel); // problem type and number of elements
	sprintf(sfnm,"%s/%s","./results",ifnm);
	
	return sptr;
}

void resoutput(double* res)
{
        int nnode = msh.nnode;
        //double **s = msh.s;
        double *s = msh.s;
	int dim = msh.dim;
        int node;
	char rfnm[50];
        FILE *resfil;
        
        strcpy(rfnm,savename());
        strcat(rfnm,".res");

	printf("Saving: %s\n",rfnm);
        resfil = fopen(rfnm,"w");
	fprintf(resfil,"%d\n",nnode);
        for (node=0;node<nnode;node++)
        {
                //fprintf(resfil,"%.15f %.15f %.15f %.15f\n",s[i][0],s[i][1],s[i][2],res[i]);
                fprintf(resfil,"%.15f %.15f %.15f %.15f\n",s[node*dim+0],s[node*dim+1],s[node*dim+2],res[node]);
        }
        fclose(resfil);
}

void erroroutput(double nerror,double* error)
{
        int nnode = msh.nnode;
        double *s = msh.s;
        int dim = msh.dim;
        int node;
        FILE *efil;
        char efnm[50];

        strcpy(efnm,savename());
        strcat(efnm,".err");

	printf("Saving: %s\n",efnm);
        efil = fopen(efnm,"w");
	fprintf(efil,"%d\n",nnode);
        fprintf(efil,"%.15f\n",nerror);
        for (node=0;node<nnode;node++)
        {
                fprintf(efil,"%.15f %.15f %.15f %.15f\n",s[node*dim+0],s[node*dim+1],s[node*dim+2],error[node]);
        }
        fclose(efil);
}
      

void vtkoutput(double* res)
{
        int nnode = msh.nnode;
        int nel = msh.nel;
        int knode = msh.knode;
	int dim = msh.dim;
        //int **icon = msh.icon;
        int *icon = msh.icon;
        //double **s = msh.s;
        double *s = msh.s;
	int *region = msh.region;
        int node,el,i;
	
        char vfnm[50];
        strcpy(vfnm,savename());
        strcat(vfnm,".vtk");
	
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

        fprintf(vtkfil,"CELL_TYPES %d\n",nel);
        for (el=0;el<nel;el++){
                fprintf(vtkfil,"%d\n",cell_type);
        }
        fprintf(vtkfil,"\n");
        fprintf(vtkfil,"POINT_DATA %d\n",nnode);
        fprintf(vtkfil,"SCALARS %s %s %d\n","resultData","float",1);
        fprintf(vtkfil,"LOOKUP_TABLE %s\n","testData");
        for (i=0;i<nnode;i++){
                fprintf(vtkfil,"%2.6f\n",res[i]);
        }
        fprintf(vtkfil,"\n");

	// print index for each part
	fprintf(vtkfil,"CELL_DATA %d\n",nel);
        fprintf(vtkfil,"SCALARS %s %s %i\n","region","int",1);
        fprintf(vtkfil,"LOOKUP_TABLE %s\n","mat");
        for (el=0;el<nel;el++)
        {
                fprintf(vtkfil,"%10d\n",region[el]);
	}
        fprintf(vtkfil,"\n");

        fclose(vtkfil);

}

