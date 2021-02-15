#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/resource.h>
#include <string.h>
#include <stdbool.h>
#include <unistd.h>


typedef struct
{
	int nnode,nel,dim,knode,neq;
	int *icon; // connectivity matrix
	double *s; // nodes
	int *bdflag; // boundary flag
	int *region; // region index
	//double *vol; // element volume - really useful to have
	int *neigh; // element neighbour indices

	// optional
	int nface; // number of internal face node
	int *iconf; // connectivity matrix for internal face
	int *bdflagf; // boundary flag for internal face
	double *lscon; // lengthscale constraint
} mesh;

// define list structure as 'list_node'
typedef struct node list_node;
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


// Define global variables here
extern mesh msh;
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
extern void listfree(list_node *);
extern void freemesh(mesh *);
extern double volume(double *,double *,double *,double *);
extern double evolume(int);
