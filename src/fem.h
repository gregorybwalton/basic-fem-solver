#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/resource.h>
#include <string.h>
#include <stdbool.h>
#include <unistd.h>

#include "tetgen.h"

// Define global variables here
struct mesh
{
	int nnode,nel,dim,knode,neq;
	int *icon;
	double *s;
	int *bdflag;
	int *region;
	double *lscon; // lengthscale constraint
};

typedef struct node list_node; // define list structure as 'list_node'
struct node
{
	int na; // neighbouring node
	double val; // stiffness matrix
	list_node *next;
};

struct sysmat
{
	int nzeros; // total number of non-zeros
	int *nnzeros; // a vector of the number of non-zeros in each row
	list_node **head; // sparse coefficient matrix linked list
	double *load; // RHS (b) vector 
	double *sol; // solution (x) vector
};

struct gauss
{
        double *w;
        double *crd;
};


extern struct mesh msh;
extern struct sysmat spstiff;
extern int soltype; // define the soltype
extern char sfnm[50]; // save file name
extern char ifnm[50]; // input file name
extern char idir[50]; // input file directory
//double *stiff, *load, *sol;
extern int nels,*icontmp;

// Define global functions here
extern double phiFunc(double,double,double);
extern double det3(double,double,double,double,double,double,double,double,double);
extern double loadFunc(double,double,double);
extern char *savename(void);
extern void printlist(list_node **,int,int);
extern void writelist(list_node **,int,int,char *);
extern void listdel(list_node **,int);
