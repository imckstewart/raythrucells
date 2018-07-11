#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <gsl/gsl_rng.h>
#include "cpgplot.h"

#include "../src/raythrucells.h"
#include "../src/meshtocube.h"
#include "../src/second_order.h"
#include "icosahedron.h"
#include "dissected_cube.h"
#include "funcs_to_test.h"

#ifndef TRUE
#define TRUE                   (_Bool)1
#endif

#ifndef FALSE
#define FALSE                  (_Bool)0
#endif

// Declarations of functions in test_prints.c:
void print_face(char *spaces, const int numDims, faceType face);

// Declarations of functions in test_utils.c:
double	flabs(double inValue);
_Bool	valuesAreClose(double valA, double valB, double epsilonMultiplier);
void	calcCrossProduct(double *vecA, double *vecB, double *vecResult);
void	matrixMultiply(const int numRowsOut, const int numColsOut\
  , const int numTermsInSum, double a[numRowsOut][numTermsInSum]\
  , double b[numTermsInSum][numColsOut], double result[numRowsOut][numColsOut]);
void	rotate(double aboutXDeg, double aboutYDeg, double aboutZDeg\
  , double *vertexCoords, unsigned long numPoints, int numDims);
void	project(const int);
double	testFunctionGauss(const double x, const double y, const double z);
double	testFunctionLinear(const double x, const double y, const double z);
void	doTestPlot(double *yValues, const int numValues);
void	freeGridPoints(struct gridPoint *gps, unsigned long numVertices);
void	getExampleCells(int *numDims, int *numElements, struct simplex **cells\
  , unsigned long *numCells, edgeType **edges, unsigned long *numEdges\
  , double **vertexCoords, double **vertexValues, struct gridPoint **gps\
  , unsigned long *numVertices);
void	getExampleCellsAndRay(int *numDims, int *numElements, struct simplex **cells\
  , unsigned long *numCells, edgeType **edges, unsigned long *numEdges\
  , double **vertexCoords, double **vertexValues, unsigned long *numVertices\
  , cellChainType *cellChain, double **rayOrigin, double **rayDir);
void	getExampleCells2ndOrder(baryBuffType *baryBuff, struct simplex **cells\
  , unsigned long *numCells, edgeType **edges, double **midEdgeValues\
  , unsigned long *numEdges, double **vertexCoords, double **vertexValues\
  , struct gridPoint **gps, unsigned long *numVertices);
void	getExampleCellsAndRay2ndOrder(baryBuffType *baryBuff, struct simplex **cells\
  , unsigned long *numCells, edgeType **edges, double **midEdgeValues\
  , unsigned long *numEdges, double **vertexCoords, double **vertexValues\
  , unsigned long *numVertices, cellChainType *cellChain\
  , double **rayOrigin, double **rayDir);

// Declarations of functions in testfunc.c:
int	check_calcEdgeVertexIndices(void);

int	check_gramSchmidt(void);
int	check_getNextEdgeSet(void);
int	check_calcDotProduct(void);
int	check_calcCellCentres(void);
int	check_gridPointsAreNeighbours(void);
int	check_edgesFormACell(void);
int	check_cellVerticesMatch(void);
int	check_addRawCell(void);
int	check_getCellsFromGrid(void);
int	check_getEdges(void);
int	check_calcBaryCoords(void);

int	check_extractFace(void);
int	check_getNewEntryFaceI(void);
int	check_calcFaceInNMinus1(void);
int	check_intersectLineWithFace(void);
int	check_followGoodChain(void);
int	check_buildRayCellChain(void);
int	check_followRayThroughCells(void);

int	check_interpolateAtFace(void);
int	check_getFaceInterpsAlongRay(void);
int	check_interpOnGridAlongRay(void);
int	check_generateVoxelIndex(void);
int	check_generateNextPixelCombo(void);
int	check_cellsToHyperCube(void);

int	check_evaluate2ndOrderShapeFns(void);
int	check_interpolate2ndOrderCell(void);
int	check_getParabolicShapeFns(void);
int	check_interpolateParabolic(void);
int	check_faceBaryToCellBary(void);
int	check_fillBaryBuffValues(void);
int	check_getInterpsAlongRay(void);
int	check_setRasterFlags(void);
int	check_interpOnGridAlongRay2ndOrder(void);


void checkRayInterp(void);
void testFits(void);

// Declarations of functions in write_fits.c:
unsigned long flattenImageIndices(const int sizes[3], const int xi, const int yi, const int zi);
void write3Dfits(char *fileName, double *cube, const int sizes[3], double xLos[3], double xHis[3]);

// Declarations of normally private functions in rt_utils:
_Bool	_getNextEdgeSet(const int numDims, const int numNeigh, _Bool *start, int neighSet[numDims]);
_Bool	_gridPointsAreNeighbours(const int numDims, struct gridPoint *gpA, struct gridPoint *gpB);
_Bool	_edgesFormACell(const int numDims, struct gridPoint *gp, int neighSet[numDims]);
_Bool	_cellVerticesMatch(const int numDims, struct simplex *cellA, struct simplex *cellB);
void	_addRawCell(const int numDims, struct simplex *candidateCell, struct simplex **cells\
  , unsigned long *maxNumCells, unsigned long *numCells);

// Declarations of normally private functions in raythrucells:
faceType	_extractFace(const int numDims, double *vertexCoords, struct simplex *dc, const unsigned long dci, const int fi);
int	_getNewEntryFaceI(const int numDims, const unsigned long dci, const struct simplex *newCell);
//facePlusBasisType _calcFaceInNMinus1_old(const int numDims, const int numVertices, faceType *face);
facePlusBasisType _calcFaceInNMinus1(const int numDims, faceType *face);
intersectType _intersectLineWithFace(const int numDims, double *x, double *dir, faceType *face\
  , const double epsilon);
_Bool	_followGoodChain(const int numDims, double *x, double *dir\
  , struct simplex *cells, unsigned long dci, int cellEntryFaceI\
  , const double epsilon, faceType **facePtrs[numDims+1]\
  , double *vertexCoords, _Bool **cellVisited, cellChainType *cellChain\
  , intersectType marginalExitIntcpts[numDims+1], int *numMarginalExits);
int	_buildRayCellChain(const int numDims, double *x, double *dir\
  , struct simplex *cells, unsigned long dci\
  , int entryFaceI, int levelI, const double epsilon\
  , faceType **facePtrs[numDims+1], double *vertexCoords\
  , _Bool **cellVisited, cellChainType *cellChain);

// Declarations of normally private functions in meshtocube:
void	_interpolateAtFace(const int numCellVertices, const int numElementsPerVertex\
  , intersectType *intercept, double *vertexValues, struct simplex *cell, double *values);
void	_getFaceInterpsAlongRay(const int numDims, const int numElementsPerVertex\
  , double *vertexValues, struct simplex *cells, cellChainType *cellChain\
  , double *faceInterpValues, double *faceDistValues);
void	_interpOnGridAlongRay(const int numDims, const int numElementsPerVertex\
  , double *vertexValues, struct simplex *cells, cellChainType *cellChain\
  , const double deltaX, const int numXi, rasterType *raster, double *rasterValues);
unsigned long	_generateVoxelIndex(const int numDims, axisType axes[N_DIMS], int *pxi);
_Bool	_generateNextPixelCombo(const int numDims, axisType axes[N_DIMS], int *pxi\
  , unsigned long *ppi);

// Declarations of normally private functions in second_order:
void	_evaluate2ndOrderShapeFns(double *barys, baryBuffType *baryBuff);
void	_interpolate2ndOrderCell(baryBuffType *baryBuff, double *barys, double *result);
void	_getParabolicShapeFns(const double x, double shapeFns[3]);
double	interpolateParabolic(double *ys, double *shapeFns);
void	_faceBaryToCellBary(baryBuffType *baryBuff, intersectType *cellFaceIntcpt, double *cellBary);
void	_fillBaryBuffValues(double *vertexValues, double *midEdgeValues\
  , struct simplex *cells, const unsigned long dci, baryBuffType *baryBuff);
void	_setRasterFlags(baryBuffType *baryBuff, const double deltaX\
  , const int numXi, double *faceDistValues, cellChainType *cellChain\
  , rasterType *raster);
void	_getInterpsAlongRay(baryBuffType *baryBuff, double *vertexValues\
  , double *midEdgeValues, struct simplex *cells, cellChainType *cellChain\
  , double *faceInterpValues, double *midInterpValues);




