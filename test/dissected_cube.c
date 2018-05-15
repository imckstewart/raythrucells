/*
We have a 3D space with a simple box-like or cuboidal shape which we want to 'tile' or dissect into tetrahedra in some sort of regular fashion. This is done in the present module first by dividing the space into L*M*N same-sized sub-cuboids, then by dissecting each of the sub-cuboids into the same pattern of 6 tetrahedra. The arrangement of the tetrahedra per cuboid can be understood as follows. Take a unit cube with its edges aligned with the coordinate axes and one corner at the origin, [0,0,0]. View this cube toward the origin along the cube diagonal that intersects the origin (therefore which also passes through the [1,1,1] corner of the cube). We dissect the cube into tetrahedra by cutting the cube along the three planes which contain that cube diagonal as well as respectively the X, Y and Z axes. Let the plane that passes through the cube diagonal as well as the X axis be named the home plane. The tetrahedra are numbered anti-clockwise from this home plane.

Abutting tetrahedra share the same vertex numbers. All six have their vertex 0 at the origin and their vertex 3 at the cube corner with coordinates [1,1,1]. Point [1,0,0] is vertex 1 for tetrahedra 0 and 5, and the remaining corners are vertices numbered 2, 1, 2 etc in sequence. Three of the tetrahedra are congruent in everything including vertex numbers, and the remaining three are congruent to the rest after mirror reflection.

Faces are numbered for their opposing vertices.
*/

#include "dissected_cube.h"

/*....................................................................*/
unsigned long
_generateVertexIndex(const int numCubes[3], const int xi, const int yi, const int zi){
  unsigned long pindex;

  pindex = (xi*(numCubes[1]+1) + yi)*(numCubes[2]+1) + zi;

  return pindex;
}

/*....................................................................*/
void
dissectedCubeVertices(const int numCubes[3], double **vertexCoords\
  , unsigned long *numPoints){
  /*
The calling routine should free *vertexCoords after use.
  */

  const int numDims=3;
  int xi,yi,zi;
  unsigned long pindex;

  *numPoints = (numCubes[0]+1)*(numCubes[1]+1)*(numCubes[2]+1);
  *vertexCoords = malloc(sizeof(**vertexCoords)*(*numPoints)*numDims);

  /* Generate the point coords.
  */
  for(xi=0;xi<numCubes[0]+1;xi++){
    for(yi=0;yi<numCubes[1]+1;yi++){
      for(zi=0;zi<numCubes[2]+1;zi++){
        pindex = _generateVertexIndex(numCubes, xi, yi, zi);
        (*vertexCoords)[pindex*numDims + 0] = (double)xi;
        (*vertexCoords)[pindex*numDims + 1] = (double)yi;
        (*vertexCoords)[pindex*numDims + 2] = (double)zi;
      }
    }
  }
}

/*....................................................................*/
int
_generateCellIndex(const int numCubes[3], const int xi, const int yi\
  , const int zi, const int ti, unsigned long *ci){

  *ci = 0; /* Just to give it some value. */

  if(xi<0 || xi>=numCubes[0] || yi<0 || yi>=numCubes[1] || zi<0 || zi>=numCubes[2])
    return 1;

  *ci = ((xi*numCubes[1] + yi)*numCubes[2] + zi)*6 + ti;

return 0;
}

/*....................................................................*/
void
generateDissectedCubeMesh(const int numCubes[3], double *vertexCoords\
  , struct simplex **cells, unsigned long *numCells){
  /*
We have a 3D space with a simple box-like or cuboidal shape which we want to 'tile' or dissect into tetrahedra in some sort of regular fashion. The present function performs this. This is done first by dividing the space into L*M*N same-sized sub-cuboids, then by dissecting each of the sub-cuboids into the same pattern of 6 tetrahedra.

Arguments:
  numCubes contains L, M, N.
  */
  const int numDims=3,numEdgesPerCell=(numDims*(numDims+1))/2;
  int xi,yi,zi,ti,vi,i,j,status;
  unsigned long pindex,ci;
  /*
Key to the XYZ offsets within each cube for the 4 vertices of each of the 6 tetrahedra per cuboid:
  */
  int voi[6][4][3] = {\
    {{0,0,0},{1,0,0},{1,1,0},{1,1,1}}\
   ,{{0,0,0},{0,1,0},{1,1,0},{1,1,1}}\
   ,{{0,0,0},{0,1,0},{0,1,1},{1,1,1}}\
   ,{{0,0,0},{0,0,1},{0,1,1},{1,1,1}}\
   ,{{0,0,0},{0,0,1},{1,0,1},{1,1,1}}\
   ,{{0,0,0},{1,0,0},{1,0,1},{1,1,1}}};

  /*
Key to the neighbours (may be in an adjoining cuboid, see cnci) abutting the 4 faces of each of the 6 tetrahedra per cube:
  */
  int cni[6][4] = {\
    {2,1,5,4}\
   ,{5,0,2,3}\
   ,{4,3,1,0}\
   ,{1,2,4,5}\
   ,{0,5,3,2}\
   ,{3,4,0,1}};

  /*
Key to the cuboid direction offsets within each cube for the neighbours abutting the 4 faces of each of the 6 tetrahedra:
  */
  int cnci[6][4][3] = {\
    {{1,0,0},{0,0,0},{0,0,0},{ 0, 0,-1}}\
   ,{{0,1,0},{0,0,0},{0,0,0},{ 0, 0,-1}}\
   ,{{0,1,0},{0,0,0},{0,0,0},{-1, 0, 0}}\
   ,{{0,0,1},{0,0,0},{0,0,0},{-1, 0, 0}}\
   ,{{0,0,1},{0,0,0},{0,0,0},{ 0,-1, 0}}\
   ,{{1,0,0},{0,0,0},{0,0,0},{ 0,-1, 0}}};

  *numCells = numCubes[0]*numCubes[1]*numCubes[2]*6;
  *cells = malloc(sizeof(**cells)*(*numCells));

  /* Generate the cells.
  */
  for(xi=0;xi<numCubes[0];xi++){
    for(yi=0;yi<numCubes[1];yi++){
      for(zi=0;zi<numCubes[2];zi++){
        for(ti=0;ti<6;ti++){
          status = _generateCellIndex(numCubes, xi, yi, zi, ti, &ci);
          if(status) {
            printf("whoops! %d %d %d %d\n", xi, yi, zi, ti);
exit(1);
          }
          for(vi=0;vi<4;vi++){
            pindex = _generateVertexIndex(numCubes, xi+voi[ti][vi][0], yi+voi[ti][vi][1], zi+voi[ti][vi][2]);
            (*cells)[ci].id = ci;
            (*cells)[ci].vertx[vi] = pindex;
          }
        }
      }
    }
  }

  /* Assign neighbours:
  */
  for(xi=0;xi<numCubes[0];xi++){
    for(yi=0;yi<numCubes[1];yi++){
      for(zi=0;zi<numCubes[2];zi++){
        for(ti=0;ti<6;ti++){
          status = _generateCellIndex(numCubes, xi, yi, zi, ti, &ci);
          if(status) {
            printf("whoops 2! %d %d %d %d\n", xi, yi, zi, ti);
exit(1);
          }
          for(vi=0;vi<4;vi++){
            status = _generateCellIndex(numCubes, xi+cnci[ti][vi][0], yi+cnci[ti][vi][1], zi+cnci[ti][vi][2], cni[ti][vi], &pindex);
            if(status)
              (*cells)[ci].neigh[vi] = NULL;
            else
              (*cells)[ci].neigh[vi] = &(*cells)[pindex];
          }
        }
      }
    }
  }

  for(i=0;i<(int)(*numCells);i++){
    for(j=0;j<numEdgesPerCell;j++)
      (*cells)[i].edges[j] = 0; /* We're not going to construct the edges here. */
  }

  /* Finally, calculate the centres of the cells: */
  calcCellCentres(numDims, *numCells, vertexCoords, *cells);
}

