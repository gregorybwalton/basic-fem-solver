#include "../fem.h"

double calcGauss(double *,int,struct gauss);
int genspstruct(list_node **,int *,int *,int,int);
struct gauss preGauss(void);
void printlist(list_node **,int,int);
void writelist(list_node **,int,int,char *);
void listadd(list_node *,int,double);

void spmatrices(void)
{
	double *s = msh.s;
	int *icon = msh.icon;
	int knode = msh.knode;
	int nel = msh.nel;
	int *bdflag = msh.bdflag;
	int *region = msh.region;
	int dim = msh.dim;
	int nt = msh.ntrue;
	int nnode = msh.nnode;
	int node,el,i,j;
	double *sloc; // Local coordinates
	double sixvol;
	double ac[knode],bc[knode],cc[knode],dc[knode];
	int vidx[6] = {1,2,3,0,1,2};
	int vsgn;
	int ii,jj;
	double kij;
	double gaussSol; // Solution to gaussian quadrature approximation
	int ibd = 0;
	int jbd = 0;
	int nzeros; // non-zeros entries
	double kappa;
	
	sloc = (double *) calloc(knode*dim,sizeof(double));
        load = (double *)calloc(msh.ntrue,sizeof(double));
        sol = (double *)calloc(msh.nnode,sizeof(double));
	list_node **arrlist = (list_node **)calloc(nt,sizeof(list_node *));

	printf("Assembling sparse matrices.\n");
	nzeros = genspstruct(arrlist,icon,bdflag,nel,knode);
	printf("Non-zero entries = %d\n",nzeros);
	spstiff.nzeros = nzeros;
	spstiff.head = arrlist;
	//printlist(arrlist,nt,0);

	// only need to declare the gauss structure first	
	struct gauss gaussvals = preGauss();

	// Assigns the Dirchlet nodes (boundary nodes?), since the solution are known here
	// In the form; [K][u]=[f]
	for (node=0;node<nnode;node++)
	{
        	if (bdflag[node]==-1)
		{
                	sol[node] = phiFunc(s[node*dim+0],s[node*dim+1],s[node*dim+2]);
	        } 
        	//printf("sol[%d]=%.2f\n",i,sol[i]);
	}

	for (el=0;el<nel;el++)
	{

		for (i=0;i<knode;i++)
		{
                	for (j=0;j<dim;j++)
			{
                        	sloc[i*dim+j] = s[icon[el*knode+i]*dim+j];
	                        //printf("sloc[%d][%d] = %.5f\n",i,j,sloc[i][j]);
			}
		}
			
		for (i=0;i<knode;i++)
		{
			// calculate the cofactors
			vsgn = (int)pow(-1,i); // introduce this sign iterator
			ac[i] = vsgn*det3(sloc[vidx[i]*dim+0],sloc[vidx[i]*dim+1],sloc[vidx[i]*dim+2],\
					sloc[vidx[i+1]*dim+0],sloc[vidx[i+1]*dim+1],sloc[vidx[i+1]*dim+2],\
					sloc[vidx[i+2]+0],sloc[vidx[i+2]+1],sloc[vidx[i+2]+2]);
			bc[i] = vsgn*det3(1.0,sloc[vidx[i]*dim+1],sloc[vidx[i]*dim+2],1.0,sloc[vidx[i+1]*dim+1],sloc[vidx[i+1]*dim+2],1.0,sloc[vidx[i+2]*dim+1],sloc[vidx[i+2]*dim+2]);
			cc[i] = -vsgn*det3(sloc[vidx[i]*dim+0],1.0,sloc[vidx[i]*dim+2],sloc[vidx[i+1]*dim+0],1.0,sloc[vidx[i+1]*dim+2],sloc[vidx[i+2]*dim+0],1.0,sloc[vidx[i+2]*dim+2]);
			dc[i] = vsgn*det3(sloc[vidx[i]*dim+0],sloc[vidx[i]*dim+1],1.0,sloc[vidx[i+1]*dim+0],sloc[vidx[i+1]*dim+1],1.0,sloc[vidx[i+2]*dim+0],sloc[vidx[i+2]*dim+1],1.0);
		}


		if (region[el] == 0)
		{
			kappa = 1.0;
		}
		else if (region[el] == 1)
		{
			kappa = 1.0e3;
		}

		// volume of element
		sixvol = 1.0*det3(sloc[3],sloc[4],sloc[5],sloc[6],sloc[7],sloc[8],sloc[9],sloc[10],sloc[11])\
				-sloc[0]*det3(1.0,sloc[4],sloc[5],1.0,sloc[7],sloc[8],1.0,sloc[10],sloc[11])\
				+sloc[1]*det3(1.0,sloc[3],sloc[5],1.0,sloc[6],sloc[8],1.0,sloc[9],sloc[11])\
				-sloc[2]*det3(1.0,sloc[3],sloc[4],1.0,sloc[6],sloc[7],1.0,sloc[9],sloc[10]);
		sixvol = fabs(sixvol);
		//printf("sixvol = %.5f\n",sixvol);

	        for (j=0;j<knode;j++)
		{
	                jj = icon[el*knode+j];
        	        jbd = bdflag[jj];

                	if (jbd>=0)
			{
				gaussSol = calcGauss(sloc,j,gaussvals);
				load[jbd] = load[jbd] + gaussSol;

				for (i=0;i<knode;i++)
				{
	                                ii = icon[el*knode+i];
        	                        ibd = bdflag[ii];
				

					// where does kappa go?	
					kij = (1.0/(6.0*sixvol))*(bc[j]*bc[i]+cc[j]*cc[i]+dc[j]*dc[i]);
					//printf("kij = %.5f\n",kij);

	                                if (ibd>=0)
					{
	                                        listadd(arrlist[jbd],ibd,kij*kappa);
					}
					else if (ibd==-1)
					{
						load[jbd] = load[jbd] - sol[ii]*kij;
					}
				}
			}
		}
	}

	//printlist(arrlist,nt);
        printf("--------------------\n");
}

int genspstruct(list_node **arrlist, int *icon, int *bdflag, int nel, int knode)
{
	// generate the sparse stiffness matrix structure
	int el,i,j;
        int nn,nm;
        int nzcount = 0; // non-zero count
	int ccount = 0; // column count
        int find;
        list_node *current;
	//list_node *new_node;

	for (el=0;el<nel;el++)
        {
                for (i=0;i<knode;i++)
                {
                        nn = bdflag[icon[el*knode+i]];
			//printf("bdflag[%d] = %d\n",icon[el*knode+i],nn);

			if (nn>=0)
			{
				// if empty
				if (arrlist[nn] == NULL)
				{
					arrlist[nn] = (list_node *)malloc(sizeof(list_node));
					arrlist[nn]->na = nn;
                                       	arrlist[nn]->val = 0.0;
                                        arrlist[nn]->next = NULL;
                                        nzcount++;
				}

	                        for (j=0;j<knode;j++)
	                        {
	                                find = 0;// used to check the occurance
					current = arrlist[nn];
	                                nm = bdflag[icon[el*knode+j]];
		
					if (nm>=0)
					{
						//printf("current->na = %d\n",current->na);
						if (current->na != nm)
	                                        {
							ccount = 0;
							while (current->next != NULL)
	                                                {
								ccount++;
								if (current->next->na == nm)
		                                                {
		                                                	find = 1;
		                                                        break;
		                                                }
		                                                current = current->next;
							}
							//printf("Row count = %d\n",ccount);

		                                        if (find == 0)
		                                        {
								current->next = (list_node *)malloc(sizeof(list_node)); 
								current->next->na = nm;
								current->next->val = 0.0;
								current->next->next = NULL;
								nzcount++;
		                                        }
	                                        }
		                        }
				}
			}
		}
	}
	return nzcount;
}


void heatscaling(scal,s,icon,nel,knode,dim)
// Generates the heat scaling for adaptivity target
double *scal,*s;
int *icon;
int nel,knode,dim;
{
	double *sloc;
	sloc = (double *) calloc(knode*dim,sizeof(double));
	int el,n,d;
	int iint;
	for (el=0;el<nel;el++)
	{
		iint = 0;
		for (n=0;n<knode;n++)
		{
			//printf("%d, coor[%d] = ",el,n);
			for (d=0;d<dim;d++)
			{
				sloc[n*dim+d] = s[icon[el*knode+n]*dim+d];
				//printf("%.5f, ",n*dim+d,sloc[n*dim+d]);
			}
			//printf("\n");
			
			if ((sloc[n*dim]>=0.25 & sloc[n*dim]<=0.75) & (sloc[n*dim+1]>=0.25 & sloc[n*dim+1]<=0.75) & (sloc[n*dim+2]>=0.25 & sloc[n*dim+2]<=0.75))
			{
				iint = 1;
			}
			else
			{
				iint = 0;
			}
		}
	
		if (iint == 1)
		{
			scal[el] = 1.0e5;
		}
		else
		{
			scal[el] = 0.0;
		}
		//printf("scal[%d] = %.2f\n",el,scal[el]);
	}

	free(sloc);
}
