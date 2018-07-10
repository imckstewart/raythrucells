#ifndef SECOND_ORDER_H
#define SECOND_ORDER_H

#include <stdio.h>
#include <math.h>

#include "dims.h" /* Should define the macro N_DIMS which specifies the number of spatial dimensions. Note that this is a _maximum_ value, the actual value 'numDims' which is passed via the function interfaces may be <= N_DIMS. */
#include "rtc_types.h"

#define SCO_ERR_TOO_MANY_VERT	1
#define SCO_ERR_BAD_CI		2

#ifndef TRUE
#define TRUE                   (_Bool)1
#endif

#ifndef FALSE
#define FALSE                  (_Bool)0
#endif


void	interpOnGridAlongRay2ndOrder(baryBuffType *baryBuff, double *vertexValues\
  , double *midEdgeValues, struct simplex *cells, const double deltaX\
  , cellChainType *cellChain, const int numXi, rasterType *raster\
  , double *rasterValues);

#endif /* SECOND_ORDER_H */

