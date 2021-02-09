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
                        fprintf(fl,"%d %d ",node,curr->na);
                        fprintf(fl,"%.5e\n",curr->val);
                        curr = curr->next;
                }
        }

        printf("Written list: %s\n",lfnm);

}

void listadd(list_node* curr,int n,double k)
// Add an scalar to an element of the linked list
// Used for adding the local stiffness element to global sparse matrix
{
        while (curr->na != n)
        {
                curr = curr->next;
        }

        curr->val += k;
}

void listdel(list_node** list,int n)
// deletes a specific item from the list
// indentify that node with n
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

void listarrayfree(list_node **list,int neq)
// free entire list array
// not tested
{
	int i;

	for (i=0;i<neq;i++)
	{
		listfree(list[i]);
	}
	//free(list); // probs need to be done outside of function
	return;
}

void listfree(list_node *head)
// free each list in the array
{
	list_node *tmp;

	while (head != NULL)
	{
		tmp = head->next;
		free(head);
		head = tmp;
	}
	return;
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

list_node** copylist(list_node **head,int nt)
// copies list to output list
{
	list_node **newhead = (list_node **)calloc(nt,sizeof(list_node *));
	list_node *newcurrent;
	list_node *current;
	int ii;

	for (ii=0;ii<nt;ii++)
	{
		newhead[ii] = (list_node *)malloc(sizeof(list_node));
		newhead[ii]->na = head[ii]->na;
		newhead[ii]->val = head[ii]->na;
		newhead[ii]->next = NULL;
		newcurrent = newhead[ii];
		current = head[ii];
		
		while (current->next != NULL)
		{
			newcurrent->next = (list_node *)malloc(sizeof(list_node));
			newcurrent->next->na = current->next->na;
			newcurrent->next->val = current->next->val;
			newcurrent->next->next = NULL;
		}
	}
	return newhead;
}
