#include "fem.h"

void writenode(char *);
void writeele(char *);
void writeface(char *);

void writetetgen(char *fnm)
{
	
	writenode(fnm);
	writeele(fnm);
	writeface(fnm);

	return;
}


void writenode(char *fnm)
{
	FILE *nfil;
	char nfnm[50];
	int n,i,abt;

	snprintf(nfnm,sizeof(nfnm),"%s%s",fnm,".node");
	nfil = fopen(nfnm,"w");
	fprintf(nfil,"%d %d %d %d\n",msh.nnode,msh.dim,0,1);

	for (n=0;n<msh.nnode;n++)
	{
		if (msh.bdflag[n] > -1)
		{
			abt = 0;
		}
		else if (msh.bdflag[n] == -1)
		{
			abt = -1;
		}
		fprintf(nfil,"%3d ",n+1);
		for (i=0;i<msh.dim;i++) fprintf(nfil,"%2.5f ",msh.s[n*msh.dim+i]);
		fprintf(nfil,"%3d\n",abt);
	}
	fclose(nfil);
	return;
}


void writeele(char *fnm)
{
	FILE *efil;
	int el,i;
	char efnm[50];

	snprintf(efnm,sizeof(efnm),"%s%s",fnm,".ele");
	efil = fopen(efnm,"w");

	fprintf(efil,"%d %d %d\n",msh.nel,msh.knode,1);
	for (el=0;el<msh.nel;el++)
	{
		fprintf(efil,"%3d ",el+1);
		for (i=0;i<msh.knode;i++) fprintf(efil,"%5d ",msh.icon[el*msh.knode+i]+1);
		fprintf(efil,"%3d\n",msh.region[el]);
	}
	fclose(efil);
	return;
}

void writeface(char *fnm)
{
	FILE *ffil;
	int n,i;
	char ffnm[50];

	snprintf(ffnm,sizeof(ffnm),"%s%s",fnm,".face");
	ffil = fopen(ffnm,"w");
	
	fprintf(ffil,"%d %d\n",msh.nface,1);
	int kf = msh.knode-1;
	for (n=0;n<msh.nface;n++)
	{
		fprintf(ffil,"%3d ",n+1);
		for (i=0;i<kf;i++) fprintf(ffil,"%5d ",msh.iconf[n*kf+i]+1);
		fprintf(ffil," %d\n",msh.bdflagf[n]);
	}
	fclose(ffil);
	return;
}
