/*
This generates a mesh of 20 connected tetrahedral cells inside the hull of a regular icosahedron.
*/

#include "icosahedron.h"
#include "rt_utils.h"

typedef struct {
  int vertices[4],shortAxisI,longAxisI,zeroAxisI;
} goldenRectangle;

/*....................................................................*/
void icosahedronVertices(double **vertexCoords, unsigned long *numPoints){
  /*
The calling routine should free *vertexCoords after use.

The vertices of an icosahedron are located at the vertices of three nested golden rectangles - i.e., a rectangle whose long vs. short sides have the golden ratio. Each of the three rectangles is aligned such that one coordinate axis is parallel to its short sides, one to its long sides, and the remaining axis is normal to the plane of the rectangle. In the struct 'goldenRectangle' the numbers of the respective axes (0 for X, 1 for Y etc) are stored in the attributes shortAxisI, longAxisI and zeroAxisI.

The rectangles are diagrammed below, showing the directions of the axes and the icosahedron vertices which are assigned to their corners:

	3--------------2   ^ X
	|              |   |
	|      0       |   |
	|              |   o------>
	1--------------0          Z

	4--------------5   ^ Y
	|              |   |
	|      1       |   |
	|              |   o------>
	6--------------7          X

	10-------------8   ^ Z
	|              |   |
	|      2       |   |
	|              |   o------>
	11-------------9          Y
  */
  int i,j,vi,vvi;
  const int numDims=3;
  const double golden = 0.5*(1.0 + sqrt(5.0));
  int shortSign[]={-1.0,-1.0,1.0,1.0},longSign[]={-1.0,1.0,-1.0,1.0};
  goldenRectangle rectangles[3];

  *numPoints = 13; /* point 12 (the last) is in the centre. */

  i = 0;
  rectangles[i].vertices[0] = 1;
  rectangles[i].vertices[1] = 0;
  rectangles[i].vertices[2] = 3;
  rectangles[i].vertices[3] = 2;
  rectangles[i].shortAxisI = 0;
  rectangles[i].longAxisI  = 2;
  rectangles[i].zeroAxisI  = 1;

  i = 1;
  rectangles[i].vertices[0] = 6;
  rectangles[i].vertices[1] = 7;
  rectangles[i].vertices[2] = 4;
  rectangles[i].vertices[3] = 5;
  rectangles[i].shortAxisI = 1;
  rectangles[i].longAxisI  = 0;
  rectangles[i].zeroAxisI  = 2;

  i = 2;
  rectangles[i].vertices[0] = 11;
  rectangles[i].vertices[1] =  9;
  rectangles[i].vertices[2] = 10;
  rectangles[i].vertices[3] =  8;
  rectangles[i].shortAxisI = 2;
  rectangles[i].longAxisI  = 1;
  rectangles[i].zeroAxisI  = 0;

  *vertexCoords = malloc(sizeof(**vertexCoords)*numDims*(*numPoints));

  /*
Here we cycle through the 3 rectangles and load in the coordinates to the relevant vertices.
  */
  for(i=0;i<3;i++){
    /* Cycle through the 4 corners of the rectangle: */
    for(j=0;j<4;j++){
      vi = rectangles[i].vertices[j];

      vvi = numDims*vi + rectangles[i].shortAxisI;
      (*vertexCoords)[vvi] = shortSign[j]*1.0;

      vvi = numDims*vi + rectangles[i].longAxisI;
      (*vertexCoords)[vvi] = longSign[j]*golden;

      vvi = numDims*vi + rectangles[i].zeroAxisI;
      (*vertexCoords)[vvi] = 0.0;
    }
  }

  /* The last point is in the centre: */
  vi = 12;
  for(i=0;i<numDims;i++){
    (*vertexCoords)[numDims*vi + i] = 0.0;
  }
}

/*....................................................................*/
void icosahedronCells(double *vertexCoords, struct simplex **cells, unsigned long *numCells){
  /*
The pointer vertexCoords should be malloc'd and filled first, preferably via icosahedronVertices().

The calling routine should free *cells after use.
  */
  int i,j;//,k,vi,vvi;
  const int numDims=3,numEdgesPerCell=(numDims*(numDims+1))/2;

  *numCells = 20;
  *cells = malloc(sizeof(**cells)*(*numCells));

  i = 0;
  (*cells)[i].id = i;
  (*cells)[i].vertx[0] =  0;
  (*cells)[i].vertx[1] =  2;
  (*cells)[i].vertx[2] = 10;
  (*cells)[i].vertx[3] = 12;
  (*cells)[i].neigh[0] = &(*cells)[ 7];
  (*cells)[i].neigh[1] = &(*cells)[ 4];
  (*cells)[i].neigh[2] = &(*cells)[ 1];
  (*cells)[i].neigh[3] = NULL;

  i = 1;
  (*cells)[i].id = i;
  (*cells)[i].vertx[0] =  0;
  (*cells)[i].vertx[1] =  8;
  (*cells)[i].vertx[2] =  2;
  (*cells)[i].vertx[3] = 12;
  (*cells)[i].neigh[0] = &(*cells)[ 9];
  (*cells)[i].neigh[1] = &(*cells)[ 0];
  (*cells)[i].neigh[2] = &(*cells)[ 2];
  (*cells)[i].neigh[3] = NULL;

  i = 2;
  (*cells)[i].id = i;
  (*cells)[i].vertx[0] =  0;
  (*cells)[i].vertx[1] =  4;
  (*cells)[i].vertx[2] =  8;
  (*cells)[i].vertx[3] = 12;
  (*cells)[i].neigh[0] = &(*cells)[15];
  (*cells)[i].neigh[1] = &(*cells)[ 1];
  (*cells)[i].neigh[2] = &(*cells)[ 3];
  (*cells)[i].neigh[3] = NULL;

  i = 3;
  (*cells)[i].id = i;
  (*cells)[i].vertx[0] =  0;
  (*cells)[i].vertx[1] =  6;
  (*cells)[i].vertx[2] =  4;
  (*cells)[i].vertx[3] = 12;
  (*cells)[i].neigh[0] = &(*cells)[13];
  (*cells)[i].neigh[1] = &(*cells)[ 2];
  (*cells)[i].neigh[2] = &(*cells)[ 4];
  (*cells)[i].neigh[3] = NULL;

  i = 4;
  (*cells)[i].id = i;
  (*cells)[i].vertx[0] =  0;
  (*cells)[i].vertx[1] = 10;
  (*cells)[i].vertx[2] =  6;
  (*cells)[i].vertx[3] = 12;
  (*cells)[i].neigh[0] = &(*cells)[ 5];
  (*cells)[i].neigh[1] = &(*cells)[ 3];
  (*cells)[i].neigh[2] = &(*cells)[ 0];
  (*cells)[i].neigh[3] = NULL;

  i = 5;
  (*cells)[i].id = i;
  (*cells)[i].vertx[0] =  6;
  (*cells)[i].vertx[1] = 10;
  (*cells)[i].vertx[2] = 11;
  (*cells)[i].vertx[3] = 12;
  (*cells)[i].neigh[0] = &(*cells)[ 6];
  (*cells)[i].neigh[1] = &(*cells)[12];
  (*cells)[i].neigh[2] = &(*cells)[ 4];
  (*cells)[i].neigh[3] = NULL;

  i = 6;
  (*cells)[i].id = i;
  (*cells)[i].vertx[0] =  7;
  (*cells)[i].vertx[1] = 11;
  (*cells)[i].vertx[2] = 10;
  (*cells)[i].vertx[3] = 12;
  (*cells)[i].neigh[0] = &(*cells)[ 5];
  (*cells)[i].neigh[1] = &(*cells)[ 7];
  (*cells)[i].neigh[2] = &(*cells)[19];
  (*cells)[i].neigh[3] = NULL;

  i = 7;
  (*cells)[i].id = i;
  (*cells)[i].vertx[0] =  2;
  (*cells)[i].vertx[1] =  7;
  (*cells)[i].vertx[2] = 10;
  (*cells)[i].vertx[3] = 12;
  (*cells)[i].neigh[0] = &(*cells)[ 6];
  (*cells)[i].neigh[1] = &(*cells)[ 0];
  (*cells)[i].neigh[2] = &(*cells)[ 8];
  (*cells)[i].neigh[3] = NULL;

  i = 8;
  (*cells)[i].id = i;
  (*cells)[i].vertx[0] =  2;
  (*cells)[i].vertx[1] =  5;
  (*cells)[i].vertx[2] =  7;
  (*cells)[i].vertx[3] = 12;
  (*cells)[i].neigh[0] = &(*cells)[18];
  (*cells)[i].neigh[1] = &(*cells)[ 7];
  (*cells)[i].neigh[2] = &(*cells)[ 9];
  (*cells)[i].neigh[3] = NULL;

  i = 9;
  (*cells)[i].id = i;
  (*cells)[i].vertx[0] =  2;
  (*cells)[i].vertx[1] =  8;
  (*cells)[i].vertx[2] =  5;
  (*cells)[i].vertx[3] = 12;
  (*cells)[i].neigh[0] = &(*cells)[16];
  (*cells)[i].neigh[1] = &(*cells)[ 8];
  (*cells)[i].neigh[2] = &(*cells)[ 1];
  (*cells)[i].neigh[3] = NULL;

  i = 10;
  (*cells)[i].id = i;
  (*cells)[i].vertx[0] =  1;
  (*cells)[i].vertx[1] =  3;
  (*cells)[i].vertx[2] =  9;
  (*cells)[i].vertx[3] = 12;
  (*cells)[i].neigh[0] = &(*cells)[17];
  (*cells)[i].neigh[1] = &(*cells)[14];
  (*cells)[i].neigh[2] = &(*cells)[11];
  (*cells)[i].neigh[3] = NULL;

  i = 11;
  (*cells)[i].id = i;
  (*cells)[i].vertx[0] =  1;
  (*cells)[i].vertx[1] = 11;
  (*cells)[i].vertx[2] =  3;
  (*cells)[i].vertx[3] = 12;
  (*cells)[i].neigh[0] = &(*cells)[19];
  (*cells)[i].neigh[1] = &(*cells)[10];
  (*cells)[i].neigh[2] = &(*cells)[12];
  (*cells)[i].neigh[3] = NULL;

  i = 12;
  (*cells)[i].id = i;
  (*cells)[i].vertx[0] =  1;
  (*cells)[i].vertx[1] =  6;
  (*cells)[i].vertx[2] = 11;
  (*cells)[i].vertx[3] = 12;
  (*cells)[i].neigh[0] = &(*cells)[ 5];
  (*cells)[i].neigh[1] = &(*cells)[11];
  (*cells)[i].neigh[2] = &(*cells)[13];
  (*cells)[i].neigh[3] = NULL;

  i = 13;
  (*cells)[i].id = i;
  (*cells)[i].vertx[0] =  1;
  (*cells)[i].vertx[1] =  4;
  (*cells)[i].vertx[2] =  6;
  (*cells)[i].vertx[3] = 12;
  (*cells)[i].neigh[0] = &(*cells)[ 3];
  (*cells)[i].neigh[1] = &(*cells)[12];
  (*cells)[i].neigh[2] = &(*cells)[14];
  (*cells)[i].neigh[3] = NULL;

  i = 14;
  (*cells)[i].id = i;
  (*cells)[i].vertx[0] =  1;
  (*cells)[i].vertx[1] =  9;
  (*cells)[i].vertx[2] =  4;
  (*cells)[i].vertx[3] = 12;
  (*cells)[i].neigh[0] = &(*cells)[15];
  (*cells)[i].neigh[1] = &(*cells)[13];
  (*cells)[i].neigh[2] = &(*cells)[10];
  (*cells)[i].neigh[3] = NULL;

  i = 15;
  (*cells)[i].id = i;
  (*cells)[i].vertx[0] =  4;
  (*cells)[i].vertx[1] =  9;
  (*cells)[i].vertx[2] =  8;
  (*cells)[i].vertx[3] = 12;
  (*cells)[i].neigh[0] = &(*cells)[16];
  (*cells)[i].neigh[1] = &(*cells)[ 2];
  (*cells)[i].neigh[2] = &(*cells)[14];
  (*cells)[i].neigh[3] = NULL;

  i = 16;
  (*cells)[i].id = i;
  (*cells)[i].vertx[0] =  5;
  (*cells)[i].vertx[1] =  8;
  (*cells)[i].vertx[2] =  9;
  (*cells)[i].vertx[3] = 12;
  (*cells)[i].neigh[0] = &(*cells)[15];
  (*cells)[i].neigh[1] = &(*cells)[17];
  (*cells)[i].neigh[2] = &(*cells)[ 9];
  (*cells)[i].neigh[3] = NULL;

  i = 17;
  (*cells)[i].id = i;
  (*cells)[i].vertx[0] =  3;
  (*cells)[i].vertx[1] =  5;
  (*cells)[i].vertx[2] =  9;
  (*cells)[i].vertx[3] = 12;
  (*cells)[i].neigh[0] = &(*cells)[16];
  (*cells)[i].neigh[1] = &(*cells)[10];
  (*cells)[i].neigh[2] = &(*cells)[18];
  (*cells)[i].neigh[3] = NULL;

  i = 18;
  (*cells)[i].id = i;
  (*cells)[i].vertx[0] =  3;
  (*cells)[i].vertx[1] =  7;
  (*cells)[i].vertx[2] =  5;
  (*cells)[i].vertx[3] = 12;
  (*cells)[i].neigh[0] = &(*cells)[ 8];
  (*cells)[i].neigh[1] = &(*cells)[17];
  (*cells)[i].neigh[2] = &(*cells)[19];
  (*cells)[i].neigh[3] = NULL;

  i = 19;
  (*cells)[i].id = i;
  (*cells)[i].vertx[0] =  3;
  (*cells)[i].vertx[1] = 11;
  (*cells)[i].vertx[2] =  7;
  (*cells)[i].vertx[3] = 12;
  (*cells)[i].neigh[0] = &(*cells)[ 6];
  (*cells)[i].neigh[1] = &(*cells)[18];
  (*cells)[i].neigh[2] = &(*cells)[11];
  (*cells)[i].neigh[3] = NULL;

  for(i=0;i<(int)(*numCells);i++){
    for(j=0;j<numEdgesPerCell;j++)
      (*cells)[i].edges[j] = 0; /* We're not going to construct the edges here. */
  }

  /* Finally, calculate the centres of the cells: */
  calcCellCentres(numDims, *numCells, vertexCoords, *cells);
}


