#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_splinalg.h>
#include <math.h>
#include <sys/resource.h>
#include <string.h>

// Define global variables here
struct mesh
{
	int nnode,nel,dim,knode,ntrue;
	int *icon;
	double *s;
	int *bdflag;
	int *region;
};

typedef struct node list_node; // define list structure as 'list_node'
struct node
{
	int na; // neighbouring node
	double val; // stiffness matrix
	list_node *next;
};

struct spmat
{
	int nzeros;
	list_node **head;
};

struct gauss
{
        double *w;
        double *crd;
};


struct mesh msh;
struct spmat spstiff;
int soltype; // define the soltype
char sfnm[50]; // save file name
char ifnm[50]; // input file name
double *stiff, *load, *sol;

// Define global functions here
double phiFunc(double,double,double);
double det3(double,double,double,double,double,double,double,double,double);
double loadFunc(double,double,double);
char *savename(void);
