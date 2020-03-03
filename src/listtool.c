#include "fem.h"

void printlist(list_node **list, int nnode, int swh)
{

        int node;
	list_node* current;
        for (node=0;node<nnode;node++)
        {
                current = list[node];

                printf("list %d = ",node);
                while (current != NULL)
                {
                        if (swh == 0)
                        {
                                printf("%d, ",current->na);
                        }
                        else if (swh == 1)
                        {
                                printf("%.5f, ",current->val);
                        }
                        current = current->next;
                }
                printf(" NULL\n");
        }
}

void writelist(list,nnz,nnode,lfnm)
list_node **list;
int nnode, nnz;
char lfnm[50];
{

        int node;
	list_node* curr;
        FILE *fl;

        fl = fopen(lfnm,"w");
        fprintf(fl,"%d\n",nnz);

        for (node=0;node<nnode;node++)
        {
                curr = list[node];

                while (curr != NULL)
                {
                        fprintf(fl,"%d %d ",node+1,curr->na+1);
                        fprintf(fl,"%.15f\n",curr->val);
                        curr = curr->next;
                }
        }

        printf("Written list: %s\n",lfnm);

}


void listdel(list,n)
list_node **list;
int n;
{
	list_node *curr, *prev;

	prev = NULL;

	for (curr = *list; curr != NULL; prev = curr, curr = curr->next)
	{
		if (curr->na == n)
		{
			if (prev == NULL)
			{
				*list = curr->next;
			}
			else
			{
				prev->next = curr->next;
			}
			free(curr);
			return;
		}
	}
}

void listtrans(spstiff)
// Transpose a linked list
// NOT TESTED
struct spmat spstiff;
{
        int nnode = msh.ntrue;
        int nt = msh.ntrue;
        list_node **list = spstiff.head;
        int node;

        list_node **tlist = (list_node **)calloc(nt,sizeof(list_node *));

        for (node=0;node<nt;node++)
        {
                tlist[node]->val = list[node]->val;
                tlist[node]->na = list[node]->na;
        }

}

void listadd(curr,n,k)
list_node *curr;
int n;
double k;
// Add an scalar to an element of the linked list
// Used for adding the local stiffness element to global sparse matrix
{
        while (curr->na != n)
        {
                curr = curr->next;
        }

        curr->val = curr->val + k;
}

double listget(curr,n)
list_node *curr;
int n;
{
	while (curr->na != n)
	{
		if (curr == NULL)
		{
			return 0.0;
		}
		curr = curr->next;
	}
	return curr->val;
}
