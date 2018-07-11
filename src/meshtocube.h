#ifndef MESHTOCUBE_H
#define MESHTOCUBE_H

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_interp.h>

#include "dims.h" /* Should define the macro N_DIMS which specifies the number of spatial dimensions. Note that this is a _maximum_ value, the actual value 'numDims' which is passed via the function interfaces may be <= N_DIMS. */
#include "rtc_types.h"

#ifndef TRUE
#define TRUE                   (_Bool)1
#endif

#ifndef FALSE
#define FALSE                  (_Bool)0
#endif


int	cellsToHyperCube(const int numDims, const int numElementsPerVertex\
  , double *vertexCoords, double *vertexValues, double *midEdgeValues\
  , struct simplex *cells, const unsigned long numCells, const double epsilon\
  , faceType **facePtrs[N_DIMS+1], axisType axes[N_DIMS], double **hypercube);

#endif /* MESHTOCUBE_H */

