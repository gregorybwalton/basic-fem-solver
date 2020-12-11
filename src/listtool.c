#include "fem.h"
// Contains additional functionality related to linked lists

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

void writelist(list_node** list,int nnz,int nnode,char* lfnm)
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


void listdel(list_node** list,int n)
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

void listtrans(struct sysmat spstiff)
// Transpose a linked list
// NOT TESTED
{
        int nt = msh.neq;
        list_node **list = spstiff.head;
        int node;

        list_node **tlist = (list_node **)calloc(nt,sizeof(list_node *));

        for (node=0;node<nt;node++)
        {
                tlist[node]->val = list[node]->val;
                tlist[node]->na = list[node]->na;
        }

}

double listget(list_node* curr,int n)
// Get value from the linked list
// NOT TESTED
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
