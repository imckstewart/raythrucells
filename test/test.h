#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include "cpgplot.h"

#include "../src/raythrucells.h"
#include "../src/meshtocube.h"
#include "icosahedron.h"
#include "dissected_cube.h"


void project(const int);
void makeHypercube(void);
unsigned long flattenImageIndices(const int sizes[3], const int xi, const int yi, const int zi);
void write3Dfits(char *fileName, double *cube, const int sizes[3], double xLos[3], double xHis[3]);
void doTestPlot(double *yValues, const int numValues);
void checkCubeIndexing(void);
void checkCubeIndexing2(void);
void checkRayInterp(void);

void testFits(void);

// Declarations of otherwise private functions in meshtocube etc:

void	_interpolateAlongRay(const int numDims, double *vertexValues, struct simplex *dc\
  , const double deltaX, intersectType entryIntcptFirstCell, unsigned long *chainOfCellIds\
  , intersectType *cellExitIntcpts, const int lenChainPtrs, rasterType *raster\
  , double *rasterValues, const int numXi, const int numElementsPerVertex);
_Bool	_generateNextPixelCombo(const int numDims, axisType axes[N_DIMS], int *pxi\
  , unsigned long *ppi);
