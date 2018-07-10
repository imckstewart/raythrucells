#ifndef RAYTHRUCELLS_H
#define RAYTHRUCELLS_H

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include "dims.h" /* Should define the macro N_DIMS which specifies the number of spatial dimensions. Note that this is a _maximum_ value, the actual value 'numDims' which is passed via the function interfaces may be <= N_DIMS. */
#include "rtc_types.h"

/* Error codes:
*/
#define RTC_ERR_SVD_FAIL	0
#define RTC_ERR_NON_SPAN	1
#define RTC_ERR_LU_DECOMP_FAIL	2
#define RTC_ERR_LU_SOLVE_FAIL	3
#define RTC_ERR_TOO_MANY_ENTRY	4
#define RTC_ERR_BUG		5
#define RTC_ERR_OLD_NOT_FOUND	6
#define RTC_ERR_BAD_CHAIN	7


#define RTC_MSG_STR_LEN		80
#define RTC_BUFFER_SIZE		1024

#ifndef TRUE
#define TRUE                   (_Bool)1
#endif

#ifndef FALSE
#define FALSE                  (_Bool)0
#endif


/* Functions in rtc_utils.c:
*/
void	rtcError(int errCode, char *message);
void	gramSchmidt(const int numDims, const int numAxes, double rawAxes[N_DIMS][N_DIMS],  double orthoAxes[N_DIMS][N_DIMS]);
double	calcDotProduct(const int numDims, double *vecA, double *vecB);
void	calcCellCentres(const int numDims, const unsigned long numCells, double *vertexCoords, struct simplex *dc);
void	getCellsFromGrid(const int numDims, struct gridPoint *gp, const unsigned long numPoints, struct simplex **cells, unsigned long *numCells);
void	getEdges(const int numDims, struct simplex *cells, const unsigned long numCells, edgeType **edges, unsigned long *numEdges);
void	calcBaryCoords(const int numDims, double vertices[N_DIMS][N_DIMS-1], double *x, double *bary);

/* Functions in raythrucells.c:
*/
int	followRayThroughCells(const int numDims, double *x, double *dir\
  , struct simplex *cells, const unsigned long numCells, const double epsilon\
  , faceType **facePtrs[numDims+1], double *vertexCoords\
  , cellChainType *cellChain);


#endif /* RAYTHRUCELLS_H */

