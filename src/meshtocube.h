#ifndef MESHTOCUBE_H
#define MESHTOCUBE_H

#include <stdio.h>
#include <math.h>

#include "dims.h" /* Should define the macro N_DIMS which specifies the number of spatial dimensions. Note that this is a _maximum_ value, the actual value 'numDims' which is passed via the function interfaces may be <= N_DIMS. */

#ifndef TRUE
#define TRUE                   (_Bool)1
#endif

#ifndef FALSE
#define FALSE                  (_Bool)0
#endif

/* In the following comments, N is short for N_DIMS, the number of spatial dimensions. */

/* This is just a convenience struct to contain a bunch of stuff which are passed very often together to functions in the 2nd-order interpolation routines. It should be initialized in initializeBaryBuf().
*/
typedef struct {
  int cellAlongRayI;
  _Bool pixelIsInCells;
} rasterType;

typedef struct {
  int numPixels;
  double origin, delta;
} axisType;


unsigned long	generateVoxelIndex(const int numDims, axisType axes[N_DIMS], int *pxi);
int	cellsToHyperCube(const int numDims, double *vertexCoords\
  , double *vertexValues, struct simplex *dc, const unsigned long numCells\
  , const double epsilon, faceType **facePtrs[N_DIMS+1], axisType axes[N_DIMS]\
  , const int numElementsPerVertex, double *midEdgeValues\
  , double **hypercube);

#endif /* MESHTOCUBE_H */

