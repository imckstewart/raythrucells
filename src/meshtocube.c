/*
This code is designed to take, as input, a connected set of simplicial cells with vertex values supplied (optionally also mid-link values?). It returns an orthonormal grid of samples, each sample being calculated as the linear (quadratic when the mid-link values are supplied) interpolation from the vertex values.

NOTE!! For the present we can only deal with convex, connected sets of cells. I.e. no ray may encounter more than 1 point of entry from non-cell to cell.

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
  , double *values, const int numElementsPerVertex){
  /*
Given a simplical cell in N dimensions, with values defined for each of the N+1 vertices, and a ray which passes through the cell, the present function uses (barycentric) linear interpolation to calculate the interpolated values at either the entry or exit point of the ray.

*** Required malloc sizes of pointers: ***
  - vertexValues: >=N*numElementsPerVertex, where N is the total number of points (the maximum value of cell->vertx for all cells).
  - values: numElementsPerVertex.
  */
  int faceVi,cellVi,ei;
  unsigned long gi,ggi;

  for(ei=0;ei<numElementsPerVertex;ei++)
    values[ei] = 0.0;

  faceVi = 0;
  for(cellVi=0;cellVi<numCellVertices;cellVi++){
    if(cellVi!=intercept.fi){
      gi = cell->vertx[cellVi];
      for(ei=0;ei<numElementsPerVertex;ei++){
        ggi = gi*numElementsPerVertex + ei;
        values[ei] += vertexValues[ggi]*intercept.bary[faceVi];
      }
      faceVi++;
    }
  }
}

/*....................................................................*/
void
_interpolateAlongRay(const int numDims, double *vertexValues, struct simplex *dc\
  , const double deltaX, intersectType entryIntcptFirstCell, unsigned long *chainOfCellIds\
  , intersectType *cellExitIntcpts, const int lenChainPtrs, rasterType *raster\
  , double *rasterValues, const int numXi, const int numElementsPerVertex){
  /*
Given the barycentric coordinates L_ of a point x_ within a simplicial cell, plus the values of some function Y at the vertices, a linear interpolation to estimate Y(x_) is given by

	        __N+1
	        \
	Y(x_) =  >    L_i * Y_i.
	        /_i=1

Here N is the number of spatial dimensions; a simplex in such a space has N+1 vertices. In the present case however we already know the barycentric coordinates (in 1 fewer dimensions) in the cell faces for the entry and exit points, thus we easily obtain the Y values for these two points, then interpolate simply along the line segment between them.

We calculate here the interpolated values for an evenly-spaced sequence of locations along the ray. The coordinate, which we label 'x', is distance along the ray from its origin. The face intercept objects supplied in the arguments entryIntcptFirstCell and cellExitIntcpts already record the value of this coordinate at each intercept point. The number of samples is supplied in argument 'numXi' and their spacing by 'deltaX'.

*** Required malloc sizes of pointers: ***
  - vertexValues: >=N*numElementsPerVertex, where N is the total number of points (the maximum value of cell->vertx for all cells).
  - dc: >= the maximum value in chainOfCellIds.
  - chainOfCellIds, cellExitIntcpts: lenChainPtrs
  - raster: numXi
  - rasterValues: numXi*numElementsPerVertex
  */
  const int numCellVertices=numDims+1;
  int xi,startXi,finisXi,si,previousSi=-1,ei,faceToggleI=0;
  double x,faceValues[2][numElementsPerVertex],fracScale,fracDist;
  unsigned long entryDci,exitDci;
  intersectType entryInt, exitInt;

  for(xi=0;xi<numXi;xi++){
    raster[xi].pixelIsInCells = FALSE; /* default - signals that the pixel is outside the cell mesh. */
    for(ei=0;ei<numElementsPerVertex;ei++)
      rasterValues[xi*numElementsPerVertex+ei] = 0.0; /* just to stop it flapping loose. */
    raster[xi].cellAlongRayI = 0; /* just to stop it flapping loose. */
  }

  if(lenChainPtrs<=0)
return;

  startXi = (int)ceil(entryIntcptFirstCell.dist/deltaX);
  if(startXi >= numXi)
return;

  if(startXi<0) startXi = 0;

  /*
Each cell (i.e. each simplex) has an ID integer. The IDs of the cells which the ray traverses (given of course in the same order the ray encounters them) are supplied in the argument 'chainOfCellIds'. We eventually need the ID of each cell we are going to interpolate within, but it is more convenient to store the index (called below 'si') of the cell in the list of cells given in 'chainOfCellIds'. Below we obtain this index for each 'xi'th raster pixel and store it in raster[xi].cellAlongRayI.
  */
  si = 0;
  for(xi=startXi;xi<numXi;xi++){
    x = deltaX*xi;

    while(si<lenChainPtrs && x>=cellExitIntcpts[si].dist)
      si++;

    if(si>=lenChainPtrs)
  break;

    raster[xi].cellAlongRayI = si;
    raster[xi].pixelIsInCells = TRUE;
  }

  if(xi<numXi)
    finisXi = xi;
  else
    finisXi = numXi-1;

  /* Now we interpolate for each pixel of the raster which falls within a cell:
  */
  for(xi=startXi;xi<=finisXi;xi++){
    si = raster[xi].cellAlongRayI; /* for brevity. */

    if(xi==startXi || si!=previousSi){
      /* In either case we need to set entry and exit quantities: namely the intercept and the cell ID it's face number relates to; the interpolated value[s] at the exit face; and (if we are at the starting pixel, or have jumped a cell between raster pixels) the interpolated value at the entry cell.
      */
      if(si==0){
        entryDci = chainOfCellIds[si]; /* points to the cell the intercept is with respect to (in this case the present cell) */
        entryInt = entryIntcptFirstCell;
      }else if(si-1==previousSi){ /* we have not skipped a cell going from one raster pixel to another */
        entryDci = exitDci;
        entryInt = exitInt;
      }else{ /* we have skipped a cell, or started deeper in the mesh than the first cell along the ray vector. */
        entryDci = chainOfCellIds[si-1]; /* points to the cell the intercept is with respect to (in this case the previous cell) */
        entryInt = cellExitIntcpts[si-1];
      }

      exitDci = chainOfCellIds[si];
      exitInt = cellExitIntcpts[si];

      if(xi==startXi || (si!=previousSi && si-1!=previousSi)){
        /* We are either at the starting pixel, or have jumped a cell between raster pixels. */
        _interpolateAtFace(numCellVertices, entryInt, vertexValues, &dc[entryDci], faceValues[faceToggleI], numElementsPerVertex);
      }

      _interpolateAtFace(numCellVertices, exitInt, vertexValues, &dc[exitDci], faceValues[1-faceToggleI], numElementsPerVertex);

      previousSi = si;
      faceToggleI = 1 - faceToggleI;
      fracScale = 1.0 / (exitInt.dist - entryInt.dist); //**** should we test somewhere that it is >0?
    }

    x = deltaX*xi;
    fracDist = (x - entryInt.dist) * fracScale;
    for(ei=0;ei<numElementsPerVertex;ei++)
      rasterValues[xi*numElementsPerVertex+ei]\
        =        fracDist *faceValues[  faceToggleI][ei]\
        + (1.0 - fracDist)*faceValues[1-faceToggleI][ei];
    /* Note that at this point faceValues[faceToggleI] should be the exit-face values and the others the entry-face values. */
  }
}

/*....................................................................*/
unsigned long
generateVoxelIndex(const int numDims, axisType axes[N_DIMS], int *pxi){
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
    *ppi = generateVoxelIndex(numDims, axes, pxi);
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

  int status=0,di,subPlanePxi[numDims-1],adi,rtcStatus=0,lenChainPtrs=0,xi,ei,i;
  unsigned long numPointsInCube,ppi,subPlanePpi;
  double rayOrigin[numDims],rayDir[numDims],*rasterValues=NULL;
  intersectType entryIntcptFirstCell,*cellExitIntcpts=NULL;
  unsigned long *chainOfCellIds=NULL;
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
      , numCells, epsilon, facePtrs, vertexCoords, &entryIntcptFirstCell, &chainOfCellIds\
      , &cellExitIntcpts, &lenChainPtrs);

    if(rtcStatus!=0){
      free(rasterValues);
      free(raster);
return rtcStatus;
    }

#ifdef WITH2ORDER
    if(midEdgeValues==NULL){
      _interpolateAlongRay(numDims, vertexValues, cells, axes[fixedAxisI].delta\
        , entryIntcptFirstCell, chainOfCellIds, cellExitIntcpts, lenChainPtrs\
        , raster, rasterValues, axes[fixedAxisI].numPixels, numElementsPerVertex);

    }else
      interpolateAlongRay2ndOrder(vertexValues, cells, axes[fixedAxisI].delta\
        , entryIntcptFirstCell, chainOfCellIds, cellExitIntcpts, lenChainPtrs\
        , raster, rasterValues, axes[fixedAxisI].numPixels, baryBuff, midEdgeValues);
#else
    _interpolateAlongRay(numDims, vertexValues, cells, axes[fixedAxisI].delta\
      , entryIntcptFirstCell, chainOfCellIds, cellExitIntcpts, lenChainPtrs\
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
  free(chainOfCellIds);
  free(cellExitIntcpts);

#ifdef WITH2ORDER
  if(midEdgeValues!=NULL)
    freeBaryBuf(baryBuff);
#endif

  return status;
}



