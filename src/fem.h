#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/resource.h>
#include <string.h>
#include <stdbool.h>
#include <unistd.h>

// Define global variables here
struct mesh
{
	int nnode,nel,dim,knode,neq;
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


struct mesh msh;
struct sysmat spstiff;
int soltype; // define the soltype
char sfnm[50]; // save file name
char ifnm[50]; // input file name
//double *stiff, *load, *sol;

// Define global functions here
double phiFunc(double,double,double);
double det3(double,double,double,double,double,double,double,double,double);
double loadFunc(double,double,double);
char *savename(void);
void printlist(list_node **,int,int);
void writelist(list_node **,int,int,char *);
void listdel(list_node **,int);
