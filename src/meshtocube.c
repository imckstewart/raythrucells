/*
This code is designed to take, as input, a connected set of simplicial cells with vertex values supplied (optionally also mid-link values?). It returns an orthonormal grid of samples, each sample being calculated as the linear (quadratic when the mid-link values are supplied) interpolation from the vertex values.

NOTE!! For the present we can only deal with a single convex, connected set of cells. I.e. no ray may encounter more than 1 point of entry from non-cell to cell.

In the comments in this module, N is short for numDims, the actual number of spatial dimensions in use, rather than the maximum number of dimensions N_DIMS, which some arrays are malloc'd to in order to obviate the need for constant mallocs and frees of the same size blocks.
*/

#include "raythrucells.h"
#include "meshtocube.h"
#ifdef WITH2ORDER
#include "second_order.h"
#endif

/*....................................................................*/
void
_interpolateAtFace(const int numCellVertices, intersectType intercept\
  , double *vertexValues, struct simplex *cell\
  , double *returnedValues, const int numElementsPerVertex){
  /*
Given a simplical cell in N dimensions, with values defined for each of the N+1 vertices, and a ray which passes through the cell, the present function uses (barycentric) linear interpolation to calculate the interpolated value(s) at either the entry or exit point of the ray.

*** Required malloc sizes of pointers: ***
  - vertexValues: >=M*numElementsPerVertex, where M is the total number of points (the maximum value of cell->vertx for all cells).
  - returnedValues: numElementsPerVertex.
  */
  int faceVi,cellVi,ei;
  unsigned long gi,ggi;

  for(ei=0;ei<numElementsPerVertex;ei++)
    returnedValues[ei] = 0.0;

  /* NOTE!! That intercept.fi must be in the range [0,numCellVertices-1]. */

  faceVi = 0;
  for(cellVi=0;cellVi<numCellVertices;cellVi++){
    if(cellVi!=intercept.fi){
      gi = cell->vertx[cellVi];
      for(ei=0;ei<numElementsPerVertex;ei++){
        ggi = gi*numElementsPerVertex + ei;
        returnedValues[ei] += vertexValues[ggi]*intercept.bary[faceVi];
      }
      faceVi++;
    }
  }
}

/*....................................................................*/
void
_getFaceInterpsAlongRay(const int numDims, double *vertexValues, struct simplex *cells\
  , cellChainType *cellChain\
  , double *faceInterpValues, double *faceDistValues, const int numElementsPerVertex){
  /*
The ray passes through M cells (where M is short for the input argument 'lenChainPtrs'). The M+1 face intercepts are already stored in the arguments 'entryIntcpt' plus the length-M list 'exitIntcpts'. The task of the present routine is twofold: (i) (for convenience) to return in 'faceDistValues' the 'dist' scalars, i.e. the distance along the ray from its origin to the point of face interception, for all M+1 face intercepts; (ii) to return in 'faceInterpValues' the M+1 values interpolated from the vertices of each intercepted face.

The pointers 'faceInterpValues' and 'faceDistValues' are expected to be malloc'd before the function is called, to the following sizes:
  - faceInterpValues: >=(M+1)*numElementsPerVertex.
  - faceDistValues: M+1.
  */

  const int numCellVertices=numDims+1;
  double *interpValuesThisFace=NULL;
  int i,ei;
  unsigned long dci;

  interpValuesThisFace = malloc(sizeof(*interpValuesThisFace)*numElementsPerVertex);

  i = 0;
  faceDistValues[i] = cellChain->entryIntcpt.dist;
  dci = cellChain->cellIds[i];
  _interpolateAtFace(numCellVertices, cellChain->entryIntcpt, vertexValues, &cells[dci], interpValuesThisFace, numElementsPerVertex);
  for(ei=0;ei<numElementsPerVertex;ei++)
    faceInterpValues[i*numElementsPerVertex+ei] = interpValuesThisFace[ei];

  for(i=1;i<=cellChain->nCellsInChain;i++){
    faceDistValues[i] = cellChain->exitIntcpts[i-1].dist;
    dci = cellChain->cellIds[i-1];
    _interpolateAtFace(numCellVertices, cellChain->exitIntcpts[i-1], vertexValues, &cells[dci], interpValuesThisFace, numElementsPerVertex);
    for(ei=0;ei<numElementsPerVertex;ei++)
      faceInterpValues[i*numElementsPerVertex+ei] = interpValuesThisFace[ei];
  }

  free(interpValuesThisFace);
}

/*....................................................................*/
void
_interpOnGridAlongRay(const int numDims, double *vertexValues, struct simplex *cells\
  , const double deltaX, cellChainType *cellChain, rasterType *raster\
  , double *rasterValues, const int numXi, const int numElementsPerVertex){
  /*
Given the cell-barycentric coordinates L_ of a point x_ within a simplicial cell, plus the values of some function Y at the vertices, a linear interpolation to estimate Y(x_) is given by

	        __N+1
	        \
	Y(x_) =  >    L_i * Y_i.
	        /_i=1

Here N is the number of spatial dimensions; a simplex in such a space has N+1 vertices. In the present case however since we already know the face-barycentric coordinates (i.e. in 1 fewer dimensions) in the cell faces for the entry and exit points, we can easily obtain the Y values for these two points, then interpolate simply along the line segment between them.

We calculate here the interpolated values for an evenly-spaced sequence of locations along the ray. The coordinate, which we label 'x', is distance along the ray from its origin. The face intercept objects supplied in the arguments entryIntcptFirstCell and cellExitIntcpts already record the value of this coordinate at each intercept point. The number of samples is supplied in argument 'numXi' and their spacing by 'deltaX'.

*** Required malloc sizes of pointers: ***
  - vertexValues: >=N*numElementsPerVertex, where N is the total number of points (the maximum value of cell->vertx for all cells).
  - cells: >= the maximum value in chainOfCellIds.
  - chainOfCellIds, cellExitIntcpts: lenChainPtrs
  - raster: numXi
  - rasterValues: numXi*numElementsPerVertex
  */

  double *faceInterpValues=NULL,*faceDistValues=NULL;
  gsl_interp_accel *acc;
  int xi,si,ei;
  double x,faceInterpValue0,faceInterpValue1,fracDist;

  faceInterpValues = malloc(sizeof(*faceInterpValues)*(cellChain->nCellsInChain+1)*numElementsPerVertex);
  faceDistValues   = malloc(sizeof(  *faceDistValues)*(cellChain->nCellsInChain+1));

  _getFaceInterpsAlongRay(numDims, vertexValues, cells, cellChain, faceInterpValues, faceDistValues, numElementsPerVertex);

  acc = gsl_interp_accel_alloc();

  /*
Each cell (i.e. each simplex) has an ID integer. The IDs of the cells which the ray traverses (given of course in the same order the ray encounters them) are supplied in the argument 'chainOfCellIds'. We eventually need the ID of each cell we are going to interpolate within, but it is more convenient to store the index (called below 'si') of the cell in the list of cells given in 'chainOfCellIds'. Below we obtain this index for each 'xi'th raster pixel and store it in raster[xi].cellAlongRayI.
  */
  for(xi=0;xi<numXi;xi++){
    x = deltaX*xi;

    if(x < faceDistValues[0] || x >= faceDistValues[cellChain->nCellsInChain]){
      raster[xi].cellAlongRayI = 0; /* just to stop it flapping loose. */
      raster[xi].pixelIsInCells = FALSE; /* default - signals that the pixel is outside the cell mesh. */
      for(ei=0;ei<numElementsPerVertex;ei++)
        rasterValues[xi*numElementsPerVertex+ei] = 0.0; /* just to stop it flapping loose. */

    }else{
      si = (int)gsl_interp_accel_find(acc, faceDistValues, (size_t)(cellChain->nCellsInChain + 1), x);
      raster[xi].cellAlongRayI = si;
      raster[xi].pixelIsInCells = TRUE;

      fracDist = (x - faceDistValues[si])/(faceDistValues[si+1] - faceDistValues[si]);

      for(ei=0;ei<numElementsPerVertex;ei++){
        faceInterpValue0 = faceInterpValues[ si   *numElementsPerVertex+ei];
        faceInterpValue1 = faceInterpValues[(si+1)*numElementsPerVertex+ei];

        rasterValues[xi*numElementsPerVertex+ei]\
          =        fracDist *faceInterpValue1\
          + (1.0 - fracDist)*faceInterpValue0;
      }
    }
  }

  gsl_interp_accel_free(acc);

  free(faceDistValues);
  free(faceInterpValues);
}

/*....................................................................*/
unsigned long
_generateVoxelIndex(const int numDims, axisType axes[N_DIMS], int *pxi){
  /*
This uses the FITS-like convention for flattening a multi-dimensional array, namely that the 0th axis is the least significant, i.e. that a change in pxi[0] is associated with the smallest change of ppi).
  */
  unsigned long ppi;
  int di;

  di = numDims-1;
  ppi = (unsigned long)pxi[di];
  for(di=numDims-2;di>=0;di--)
    ppi = (unsigned long)axes[di].numPixels*ppi + (unsigned long)pxi[di];

return ppi;
}

/*....................................................................*/
_Bool
_generateNextPixelCombo(const int numDims, axisType axes[N_DIMS], int *pxi\
  , unsigned long *ppi){
  /*
We have a hypercube of voxels of numDims dimensions, and want to cycle through all the combinations of pixel coordinates. Given a set of axis dimensions supplied in 'axes' and a pointer to the set of pixel coordinates 'pxi', the present function calculates and returns in 'pxi' the next combination. If it has run through all the combinations, or if the entry combination already exceeds any of the dimensions, it returns isFinished = TRUE.

Since multidimensional arrays are often most conveniently passed around in C as single-dimensional arrays, this function also calculates and returns in 'ppi' an appropriate single-dimensional index (using the FITS-like convention here that the 0th axis is the least significant, i.e. that a change in pxi[0] is associated with the smallest change of ppi).
  */
  int di;
  _Bool isFinished=FALSE;

  di = 0;
  while(di<numDims){
    pxi[di]++;

    if(pxi[di]<axes[di].numPixels)
  break;

    pxi[di] = 0;
    di++;
  }

  *ppi = 0;
  if(di>=numDims) /* This happens if we didn't break out of the earlier loop. */
    isFinished = TRUE;
  else{
    *ppi = _generateVoxelIndex(numDims, axes, pxi);
  }

return isFinished;
}

/*....................................................................*/
int
cellsToHyperCube(const int numDims, double *vertexCoords\
  , double *vertexValues, struct simplex *cells, const unsigned long numCells\
  , const double epsilon, faceType **facePtrs[N_DIMS+1], axisType axes[N_DIMS]\
  , const int numElementsPerVertex, double *midEdgeValues\
  , double **hypercube){
  /*
This function takes as input a connected mesh of simplicial cells defined by the arguments 'cells' and 'vertexValues' and returns a hypercube (i.e. an orthogonal array in the 'numDims' number of spatial dimensions) of interpolated samples. The mesh vertices, and therefore the returned samples, may be either scalar- or vector-valued, as indicated by the input argument 'numElementsPerVertex'.

The samples are taken at equally-spaced points along a set of N-1 rays which are aligned with the zeroth spatial axis.

The calling routine should free *hypercube after use.

Below is a code snippet showing an example of how to call this function:

  const int numDims=3,numElementsPerVertex=1;
  double *vertexCoords=NULL,*hypercube=NULL,*vertexValues=NULL;
  unsigned long numPoints=0,numCells=0;
  struct simplex *cells=NULL;
  const double epsilon=1.0e-6;
  axisType axes[numDims];
  int status=0;

  // Get numPoints.
  vertexCoords = malloc(sizeof(*vertexCoords)*numPoints*numDims);
  // Load vertexCoords.

  vertexValues = malloc(sizeof(*vertexValues)*numPoints*numElementsPerVertex);
  // Load vertexValues.

  // Get numCells.
  cells = malloc(sizeof(*cells)*numCells);
  // Load cells.

  // Load axes.

  status = cellsToHyperCube(numDims, vertexCoords, vertexValues, cells, numCells\
    , epsilon, NULL, axes, numElementsPerVertex, NULL, &hypercube);

  free(hypercube);
  free(vertexValues);
  free(cells);
  free(vertexCoords);



*** Required malloc sizes of input pointers: ***
  - midEdgeValues: >=M*numElementsPerVertex, where M is the total number of edges (the maximum value of cell->edges for all cells).
  */

  int status=0,di,subPlanePxi[numDims-1],rtcStatus=0,xi,ei;
  unsigned long numPointsInCube,ppi,subPlanePpi;
  double rayOrigin[numDims],rayDir[numDims],*rasterValues=NULL;
  cellChainType cellChain;
  rasterType *raster=NULL;
  _Bool finished;
#ifdef WITH2ORDER
  baryBuffType *baryBuff=NULL;
#endif
  axisType subPlaneAxes[N_DIMS];
  const int fixedAxisI=0;

  for(di=1;di<numDims;di++)
    subPlaneAxes[di-1] = axes[di];

  raster       = malloc(sizeof(*raster)*axes[fixedAxisI].numPixels);
  rasterValues = malloc(sizeof(*raster)*axes[fixedAxisI].numPixels*numElementsPerVertex);

  numPointsInCube = 1;
  for(di=0;di<numDims;di++)
    numPointsInCube *= (unsigned long)axes[di].numPixels;
//*** test for overflow!!

  *hypercube = malloc(sizeof(**hypercube)*numPointsInCube*numElementsPerVertex);

  /*
We will generate a hyperplanar set of points, i.e. a hypercube of rank 1 less than the one we are trying to fill in the present routine. This hyperplane will have 1 point at each pixel combination of the last N-1 axes of the hypercube, holding the zeroth axis pixel at 0. Each point will be the starting point for a ray, which will be followed through the set of cells by incrementing the zeroth axis pixel, with a value being calculated for each sample.
  */

  rayDir[0] = 1.0;
  for(di=1;di<numDims;di++)
    rayDir[di] = 0.0;

  rayOrigin[fixedAxisI] = axes[fixedAxisI].origin; /* for all the rays. */

  for(di=0;di<numDims-1;di++)
    subPlanePxi[di] = 0;

#ifdef WITH2ORDER
  if(midEdgeValues!=NULL)
    initializeBaryBuf(numDims, numElementsPerVertex, baryBuff);
#endif

//**** fill vertexValues if facePtrs==NULL?

  /* This thing generates all combinations of pixel values for the last N-1 axes.
  */
  ppi = 0;
  finished = FALSE;
  while(!finished){
    /* Sample at each pixel along the ray starting at the given coords for axes 1 to N-1 and extending in the direction of the 0th axis.
    */
    for(di=1;di<numDims;di++)
      rayOrigin[di] = axes[di].origin + subPlanePxi[di-1]*axes[di].delta;

    rtcStatus = followRayThroughCells(numDims, rayOrigin, rayDir, cells\
      , numCells, epsilon, facePtrs, vertexCoords, &cellChain);

    if(rtcStatus!=0){
      free(rasterValues);
      free(raster);
return rtcStatus;
    }

#ifdef WITH2ORDER
    if(midEdgeValues==NULL){
      _interpOnGridAlongRay(numDims, vertexValues, cells, axes[fixedAxisI].delta\
        , &cellChain\
        , raster, rasterValues, axes[fixedAxisI].numPixels, numElementsPerVertex);


    }else
      interpOnGridAlongRay2ndOrder(vertexValues, cells, axes[fixedAxisI].delta\
        , &cellChain\
        , raster, rasterValues, axes[fixedAxisI].numPixels, baryBuff, midEdgeValues);
#else
    _interpOnGridAlongRay(numDims, vertexValues, cells, axes[fixedAxisI].delta\
      , &cellChain\
      , raster, rasterValues, axes[fixedAxisI].numPixels, numElementsPerVertex);
#endif

    for(xi=0;xi<axes[fixedAxisI].numPixels;xi++){
      for(ei=0;ei<numElementsPerVertex;ei++){
        (*hypercube)[(ppi+xi)*numElementsPerVertex+ei] = rasterValues[xi*numElementsPerVertex+ei];
      }
    }

    finished = _generateNextPixelCombo(numDims-1, subPlaneAxes, subPlanePxi, &subPlanePpi);
    ppi = (unsigned long)axes[fixedAxisI].numPixels*subPlanePpi;
  }

  free(rasterValues);
  free(raster);
  free_cellChain(&cellChain);

#ifdef WITH2ORDER
  if(midEdgeValues!=NULL)
    freeBaryBuf(baryBuff);
#endif

  return status;
}



