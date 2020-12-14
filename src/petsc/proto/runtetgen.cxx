#include "../../fem.h"
//#include "tetgen.h"

void runtetgen()
{
	printf("\n\nRunning: runtetgen\n");
	int n,d,el,k;
	tetgenio in,out;
	int dim = msh.dim;
	int knode = msh.knode;
	//char switches[50] = "rkqaO0"; // switches being used in refinement
	tetgenbehavior switches;
	switches.refine = 1; // refines
	switches.vtkview = 1; // outputs vtks
	switches.quality = 1;
	switches.varvolume = 1;
	switches.opt_scheme = 0;
	//switches.quiet =  1; // supresses terminal output, except errors

	in.firstnumber = 0;
	in.mesh_dim = msh.dim;
	
	in.numberofpoints = msh.nnode;
	in.pointlist = new REAL[msh.nnode * msh.dim];
	in.numberofpointattributes = 1;
	in.pointattributelist = new REAL[msh.nnode];
	
        for (n=0;n<msh.nnode;n++)
        {
                for (d=0;d<dim;d++)
                {
			in.pointlist[dim*n+d] = (REAL)msh.s[dim*n+d];
		}
		if (msh.bdflag[n] > -1)
                {
			in.pointattributelist[n] = 0.0;
		}
		else
		{
			in.pointattributelist[n] = -1.0;
		}

	}

	in.numberoftetrahedra = msh.nel;
	in.tetrahedronlist = new int[msh.nel*msh.knode];
	in.tetrahedronvolumelist = new REAL[msh.nel];
	in.numberoftetrahedronattributes = 0;
	if (msh.region != NULL)
	{
		in.numberoftetrahedronattributes = 1;
		in.tetrahedronattributelist = new REAL[msh.nel];
	}

        for (el=0;el<msh.nel;el++)
        {
		in.tetrahedronvolumelist[el] = (REAL)msh.lscon[el];
                for (k=0;k<knode;k++)
                {
                        in.tetrahedronlist[knode*el+k] = msh.icon[knode*el+k];
                }
		if (in.numberoftetrahedronattributes == 1) in.tetrahedronattributelist[el] = (REAL)msh.region[el];
        }

	// addin, -i, additional vertices
	// bgmin, -m, constraint functions
	tetrahedralize(&switches,&in,&out);
	return;
}

void runtetgenexample()
{
  tetgenio in, out;
  tetgenio::facet *f;
  tetgenio::polygon *p;
  int i;

	printf("\n\nRunning: runtetgenexample\n");

in.firstnumber = 1;

 in.numberofpoints = 8;
  in.pointlist = new REAL[in.numberofpoints * 3];
  in.pointlist[0]  = 0;  // node 1.
  in.pointlist[1]  = 0;
  in.pointlist[2]  = 0;
  in.pointlist[3]  = 2;  // node 2.
  in.pointlist[4]  = 0;
  in.pointlist[5]  = 0;
  in.pointlist[6]  = 2;  // node 3.
  in.pointlist[7]  = 2;
  in.pointlist[8]  = 0;
  in.pointlist[9]  = 0;  // node 4.
  in.pointlist[10] = 2;
  in.pointlist[11] = 0;
  for (i = 4; i < 8; i++) {
    in.pointlist[i * 3]     = in.pointlist[(i - 4) * 3];
    in.pointlist[i * 3 + 1] = in.pointlist[(i - 4) * 3 + 1];
    in.pointlist[i * 3 + 2] = 12;
	}

  in.numberoffacets = 6;
  in.facetlist = new tetgenio::facet[in.numberoffacets];
  in.facetmarkerlist = new int[in.numberoffacets];

  f = &in.facetlist[0];
  f->numberofpolygons = 1;
  f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
  f->numberofholes = 0;
  f->holelist = NULL;
  p = &f->polygonlist[0];
  p->numberofvertices = 4;
  p->vertexlist = new int[p->numberofvertices];
  p->vertexlist[0] = 1;
  p->vertexlist[1] = 2;
  p->vertexlist[2] = 3;
  p->vertexlist[3] = 4;

  f = &in.facetlist[1];
  f->numberofpolygons = 1;
  f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
  f->numberofholes = 0;
  f->holelist = NULL;
  p = &f->polygonlist[0];
  p->numberofvertices = 4;
  p->vertexlist = new int[p->numberofvertices];
  p->vertexlist[0] = 5;
  p->vertexlist[1] = 6;
  p->vertexlist[2] = 7;
  p->vertexlist[3] = 8;

  f = &in.facetlist[2];
  f->numberofpolygons = 1;
  f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
  f->numberofholes = 0;
  f->holelist = NULL;
  p = &f->polygonlist[0];
  p->numberofvertices = 4;
  p->vertexlist = new int[p->numberofvertices];
  p->vertexlist[0] = 1;
  p->vertexlist[1] = 5;
  p->vertexlist[2] = 6;
  p->vertexlist[3] = 2;

  f = &in.facetlist[3];
  f->numberofpolygons = 1;
  f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
  f->numberofholes = 0;
  f->holelist = NULL;
  p = &f->polygonlist[0];
  p->numberofvertices = 4;
  p->vertexlist = new int[p->numberofvertices];
  p->vertexlist[0] = 2;
  p->vertexlist[1] = 6;
  p->vertexlist[2] = 7;
  p->vertexlist[3] = 3;

  f = &in.facetlist[4];
  f->numberofpolygons = 1;
  f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
  f->numberofholes = 0;
  f->holelist = NULL;
  p = &f->polygonlist[0];
  p->numberofvertices = 4;
  p->vertexlist = new int[p->numberofvertices];
  p->vertexlist[0] = 3;
  p->vertexlist[1] = 7;
  p->vertexlist[2] = 8;
  p->vertexlist[3] = 4;


  f = &in.facetlist[5];
  f->numberofpolygons = 1;
  f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
  f->numberofholes = 0;
  f->holelist = NULL;
  p = &f->polygonlist[0];
  p->numberofvertices = 4;
  p->vertexlist = new int[p->numberofvertices];
  p->vertexlist[0] = 4;
  p->vertexlist[1] = 8;
  p->vertexlist[2] = 5;
  p->vertexlist[3] = 1;

  in.facetmarkerlist[0] = -1;
  in.facetmarkerlist[1] = -2;
  in.facetmarkerlist[2] = 0;
  in.facetmarkerlist[3] = 0;
  in.facetmarkerlist[4] = 0;
  in.facetmarkerlist[5] = 0;

  //in.save_nodes("barin");
  //in.save_poly("barin");


        tetgenbehavior switches;
        switches.plc = 1;
        switches.minratio = 1.414;
        switches.maxvolume = 0.1;

  tetrahedralize(&switches, &in, &out);
	//tetrahedralize("pq1.414a0.1", &in, &out);

  //out.save_nodes("barout");
  //out.save_elements("barout");
  //out.save_faces("barout");

  return;
}
