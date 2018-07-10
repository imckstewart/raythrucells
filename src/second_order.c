/*
This contains the necessary routines for meshgrid to perform a 2nd-order interpolation.

In the comments in this module, N is short for numDims, the actual number of spatial dimensions in use, rather than the maximum number of dimensions N_DIMS, which some arrays are malloc'd to in order to obviate the need for constant mallocs and frees of the same size blocks.

General barycentric interpolation is essentially the same technique as the use of shape functions in Finite Element analysis. The idea is that one starts with a cell containing N fixed points (each of which may be either on the boundary of the cell or inside it) at which the values V_i are known. If one defines N shape functions Q_i within the cell, the ith shape function having the property that it is zero-valued at each of the N points except the ith, and has value unity there, then the interpolation value at point r_ is given by

	       _N
	       \
	f(r_) =  >    V_i*Q_i(r_).
	       /_i=1

For a linear interpolation, the vertices of the cell form the set of points, and each shape function Q_i(r_) is planar, and in fact equal to the ith barycentric coordinate B_i of r_. A 2nd-order interpolation requires additional V_i values at the midpoints of the edges between vertices. The shape functions here fall into two sorts, depending on whether the unity point is a cell vertex or a half-edge point. Those Q_i which have value unity at the ith vertex and zero at all the other vertices, and also at the half-edge points, are given by

	Q_i(r_) = B_i(2*B_i - 1).

For the point ij half-way between vertices i and j the appropriate shape functions, of the second sort, are

	Q_ij(r_) = 4*B_i*B_j.

The function meshtocube.c:cellsToHyperCube() takes a set of connected simplices, with values provided for the vertices, and interpolates these values onto a hypercubical grid of samples. The routine performs this interpolation by defining a marginal grid in N-1 dimensions, then following a ray from each of these grid locations in the direction of the final (the Nth, starting the count at 1) coordinate direction. The function that performs this ray sampling in meshgrid.c is called _interpOnGridAlongRay(). For 2nd-order interpolation, that function is replaced by the one called interpOnGridAlongRay2ndOrder() in the present module.
*/


#ifdef TEST
#include "../test/test.h"
#endif

#include "raythrucells.h"
#include "second_order.h"

/*....................................................................*/
void
_evaluate2ndOrderShapeFns(double *barys, baryBuffType *baryBuff){
  /*
This evaluates the 2nd-order shape functions Q_i(r_) at r_ as described in the header, as used in the function _interpolate2ndOrderCell(). The net interpolated value will then be given by

  sum(vertexValues*baryBuff->vertxShapeValues) + sum(midEdgeValues*baryBuff->edgeShapeValues)

Input arguments:
================
	- barys: this contains cell-specific barycentric coordinates of the point we want to interpolate at.

Input/output arguments:
=======================
	- baryBuff: must have been initialized via init_baryBuff().
  */
  int vi,ei;

  for(vi=0;vi<baryBuff->numVertices;vi++)
    baryBuff->vertxShapeValues[vi] = barys[vi]*(2.0*barys[vi] - 1.0);

  for(ei=0;ei<baryBuff->numEdges;ei++)
    baryBuff->edgeShapeValues[ei]\
       = 4.0*barys[baryBuff->edgeVertexIndices[ei][0]]\
            *barys[baryBuff->edgeVertexIndices[ei][1]];
}

/*....................................................................*/
void
_interpolate2ndOrderCell(baryBuffType *baryBuff, double *barys, double *result){
  /*
This evaluates the sum

	        _N
	        \
	f(r_) =  >    V_i*Q_i(r_).
	        /_i=1

as described in the module header. (Note that f and V_i are assumed here to be vector quantities, of size baryBuff->numElements.)

Input arguments:
================
	- baryBuff: must have been filled via _fillBaryBuffValues().

	- barys: this contains cell-specific barycentric coordinates of the point we want to interpolate at.

Output arguments:
=================
	- result: should be malloc'd to baryBuff->numElements.
  */
  int ei,vi,li;

  _evaluate2ndOrderShapeFns(barys, baryBuff);

  for(ei=0;ei<baryBuff->numElements;ei++)
    result[ei] = 0.0;

  for(vi=0;vi<baryBuff->numVertices;vi++)
    for(ei=0;ei<baryBuff->numElements;ei++)
      result[ei] += baryBuff->vertexValues[vi*baryBuff->numElements+ei]*baryBuff->vertxShapeValues[vi];

  for(li=0;li<baryBuff->numEdges;li++){
    for(ei=0;ei<baryBuff->numElements;ei++)
      result[ei] += baryBuff->edgeValues[li*baryBuff->numElements+ei]*baryBuff->edgeShapeValues[li];
  }
}

/*....................................................................*/
void
_getParabolicShapeFns(const double x, double shapeFns[3]){
  /*
This calculates the shape functions for interpolateParabolic().
  */

  shapeFns[0] = (x - 1.0)*(2.0*x - 1.0);
  shapeFns[1] = 4.0*x*(1.0 - x);
  shapeFns[2] = x*(2.0*x - 1.0);
}

/*....................................................................*/
double
interpolateParabolic(double *ys, double *shapeFns){
  /*
This function is supplied with values of y at the beginning, end and midpoint of the line segment; what it does is perform a 2nd-order (i.e. quadratic) interpolation to obtain an estimate of y at the given fractional displacement along the line segment.

Given y0, y1 and y2 at x0, x1 and x2, the Lagrange interpolating polynomial is

	         x-x1    x-x2         x-x0    x-x2         x-x0    x-x1
	y ~ y0*-------*------- + y1*-------*------- + y2*-------*-------.
	        x0-x1   x0-x2        x1-x0   x1-x2        x2-x0   x2-x1

If x is the fractional distance along the line segment, and [x0,x1,x2] = [0.0,0.5,1.0], then for the y values we have in hand, this reduces to

	y ~ y0*(x-1)*(2x-1) + y1*4x*(1-x) + y2*x*(2x-1).

The parts of the terms which are functions of x are called shape functions.
  */
  double result;

  result = ys[0]*shapeFns[0] + ys[1]*shapeFns[1] + ys[2]*shapeFns[2];

return result;
}

/*....................................................................*/
void
_faceBaryToCellBary(baryBuffType *baryBuff, intersectType *cellFaceIntcpt\
  , double *cellBary){
  /*
Here we convert the face-specific barycentric coordinates of the point where the ray crosses one of the faces of the cell to baycentric coordinates specific to the cell which the intercept is tied to. This consists simply of increasing the size of the bary vector by 1 element, and setting the element corresponding to the exit face to zero. The present function may be called for all exit-face intercepts and for the entry intercept for the first cell encountered.

Input arguments:
================
	- baryBuff: must have been initialized via init_baryBuff(). After that, it contains some useful integers giving the number of spatial dimensions, number of cell vertices etc.

	- cellFaceIntcpt: this contains the face-specific barycentric coordinates of the intercept between the ray and that face. The index of that face in the cell is also stored.

Output arguments:
=================
	- cellBary: should be malloc'd to baryBuff->numVertices.
  */
  int vi,vvi;

  vvi = 0;
  for(vi=0;vi<baryBuff->numVertices;vi++){
    if(vi==cellFaceIntcpt->fi){
      cellBary[vi] = 0.0;
    }else{
      cellBary[vi] = cellFaceIntcpt->bary[vvi];
      vvi++;
    }
  }
}

/*....................................................................*/
void
_fillBaryBuffValues(double *vertexValues, double *midEdgeValues\
  , struct simplex *cells, const unsigned long dci, baryBuffType *baryBuff){
  /*
'baryBuff' has been designed as a convenient structure to store a lot of things pertaining to a cell, which avoids lengthy argument lists in fucntion calls. Here the values of the interpolation function at both vertices and mid-edges of the cell cells[dci] are copied to baryBuff.

Input arguments:
================
	- vertexValues: contains samples at all cell vertices of the function to be interpolated.

	- midEdgeValues: contains samples at the mid-point of all cell edges (i.e. links between vertices) of the function to be interpolated. Note that such mid-edge values are necessary for 2nd-order interpolation.

	- cells: contains geometrical information about all the cells in the set. Mostly it is used to connect a given cell with the values in 'vertexValues' and 'midEdgeValues' which correspond to that cell's vertices and edges.

	- dci: designates the cell whose values we want to copy.

Input/output arguments:
=======================
	- baryBuff: must have been initialized via init_baryBuff(). After that, it contains some useful integers giving the number of spatial dimensions, number of cell vertices etc.


  */
  int i,ei,bvi,bli;
  unsigned long vi,li;

  for(i=0;i<baryBuff->numVertices;i++){
    vi = cells[dci].vertx[i];
    for(ei=0;ei<baryBuff->numElements;ei++){
      bvi = i*baryBuff->numElements + ei;
      baryBuff->vertexValues[bvi] = vertexValues[vi*baryBuff->numElements + ei];
    }
  }

  for(i=0;i<baryBuff->numEdges;i++){
    li = cells[dci].edges[i];
    for(ei=0;ei<baryBuff->numElements;ei++){
      bli = i*baryBuff->numElements + ei;
      baryBuff->edgeValues[bli] = midEdgeValues[li*baryBuff->numElements + ei];
    }
  }
}

/*....................................................................*/
void
_getInterpsAlongRay(baryBuffType *baryBuff, double *vertexValues\
  , double *midEdgeValues, struct simplex *cells, cellChainType *cellChain\
  , double *faceInterpValues, double *midInterpValues){
  /*
We have a ray which passes through a connected set of simplicial cells. The IDs of the cells which the ray passes through have already been worked out at this point, all the relevant info (including all the cell-face intercept locations) being stored in the argument 'cellChain'. In the present routine we want to calculated (second-order) interpolated values at each face intercept (which are returned in 'faceInterpValues') as well as the interpolated values half-way between each pair of face intercept locations (returned in 'midInterpValues').

Input arguments:
================
	- baryBuff: must have been initialized via init_baryBuff(). After that, it contains some useful integers giving the number of spatial dimensions, number of cell vertices etc.

	- vertexValues: contains samples at all cell vertices of the function to be interpolated.

	- midEdgeValues: contains samples at the mid-point of all cell edges (i.e. links between vertices) of the function to be interpolated. Note that such mid-edge values are necessary for 2nd-order interpolation.

	- cells: contains geometrical information about all the cells in the set. Mostly it is used to connect a given cell with the values in 'vertexValues' and 'midEdgeValues' which correspond to that cell's vertices and edges.

	- cellChain: contains information (chiefly cell IDs and cell-face ray intercepts) about all the cells traversed by the ray.

Output arguments:
=================
	- faceInterpValues: must be malloc'd to size (cellChain->nCellsInChain+1)*baryBuff->numElements.

	- midInterpValues: must be malloc'd to size cellChain->nCellsInChain*baryBuff->numElements.
  */

  _Bool *doFace=NULL;
  int fi,ci,ei,vi;
  double (*faceBarys)[N_DIMS+1]=NULL;
  double *interpValuesThisFace=NULL,barys[N_DIMS+1];
  unsigned long dci;

  doFace = malloc(sizeof(   *doFace)*(cellChain->nCellsInChain+1));

  for(fi=0;fi<cellChain->nCellsInChain+1;fi++)
    doFace[fi] = FALSE; /* default */

  for(ci=0;ci<cellChain->nCellsInChain;ci++){
    if(cellChain->flags[ci]){
      doFace[ci  ] = TRUE;
      doFace[ci+1] = TRUE;
    }
  }

  faceBarys = malloc(sizeof(*faceBarys)*(cellChain->nCellsInChain+1));

  fi = 0;
  if(doFace[fi])
    _faceBaryToCellBary(baryBuff, &cellChain->entryIntcpt, faceBarys[fi]);

  for(fi=1;fi<cellChain->nCellsInChain+1;fi++){
    if(doFace[fi])
      _faceBaryToCellBary(baryBuff, &cellChain->exitIntcpts[fi-1], faceBarys[fi]);
  }

  interpValuesThisFace = malloc(sizeof(*interpValuesThisFace)*baryBuff->numElements); //*** better to pass in a buffer pre-malloc'd. It would save time.

  for(fi=0;fi<cellChain->nCellsInChain+1;fi++){
    if(doFace[fi]){
      /* fill in the vertexValues and edgeValues of baryBuff.
      */
      if(fi==0)
        dci = cellChain->cellIds[fi];
      else
        dci = cellChain->cellIds[fi-1];

      _fillBaryBuffValues(vertexValues, midEdgeValues, cells, dci, baryBuff);

      _interpolate2ndOrderCell(baryBuff, faceBarys[fi], interpValuesThisFace);

      for(ei=0;ei<baryBuff->numElements;ei++)
        faceInterpValues[fi*baryBuff->numElements+ei] = interpValuesThisFace[ei];

    }else{
      for(ei=0;ei<baryBuff->numElements;ei++)
        faceInterpValues[fi*baryBuff->numElements+ei] = 0.0; /* just to stop it flapping loose. */
    }
  }

  for(ci=0;ci<cellChain->nCellsInChain;ci++){
    if(cellChain->flags[ci]){
      /* fill in the vertexValues and edgeValues of baryBuff.
      */
      _fillBaryBuffValues(vertexValues, midEdgeValues, cells, cellChain->cellIds[ci], baryBuff);

      /* interpolate the mid-interval bary coords.
      */
      for(vi=0;vi<baryBuff->numVertices;vi++)
        barys[vi] = 0.5*(faceBarys[ci][vi] + faceBarys[ci+1][vi]);

      _interpolate2ndOrderCell(baryBuff, barys, interpValuesThisFace);

      for(ei=0;ei<baryBuff->numElements;ei++)
        midInterpValues[ci*baryBuff->numElements+ei] = interpValuesThisFace[ei];

    }else{
      for(ei=0;ei<baryBuff->numElements;ei++)
        midInterpValues[ci*baryBuff->numElements+ei] = 0.0; /* just to stop it flapping loose. */
    }
  }

  free(interpValuesThisFace);
  free(faceBarys);
  free(doFace);
}

/*....................................................................*/
void
_setRasterFlags(baryBuffType *baryBuff, const double deltaX\
  , const int numXi, double *faceDistValues, cellChainType *cellChain\
  , rasterType *raster){
  /*
We have a ray which passes through a connected set of simplicial cells. We plan to calculate here the interpolated values for an evenly-spaced sequence of locations along the ray. The coordinate, which we label 'x', is distance along the ray from its origin. Here we supply in 'cellChain' a list of the cells which the ray traverses. What the present function does is check and flag which cell each of the 'numXi' sample points falls within.

Note that the present function also sets cellChain->flags for each cell that contains a raster sample point.

Input arguments:
================
	- baryBuff: must have been initialized via init_baryBuff(). After that, it contains some useful integers giving the number of spatial dimensions, number of cell vertices etc.

	- deltaX: the spacing of the samples.

	- numXi: required number of samples along the ray.

	- faceDistValues: a list of the distances along the ray at which it encounters each cell face. If there are N cells in cellChain, there will be N+1 entries in this list.

Input/output arguments:
=======================

	- cellChain: contains information (chiefly cell IDs and cell-face ray intercepts) about all the cells traversed by the ray.

Output arguments:
=================
	- raster: must have been malloc'd to size numXi. Contains the index (in the lists cellChain->cellIds etc) of the cell each sample falls within, plus a flag to indicate whether it falls within one at all.
  */

  int i,xi,si;
  gsl_interp_accel *acc;
  double x;

  for(i=0;i<cellChain->nCellsInChain;i++)
    cellChain->flags[i] = FALSE; /* default */

  acc = gsl_interp_accel_alloc();

  for(xi=0;xi<numXi;xi++){
    x = deltaX*xi;

    if(x < faceDistValues[0] || x >= faceDistValues[cellChain->nCellsInChain]){
      raster[xi].cellAlongRayI = 0; /* just to stop it flapping loose. */
      raster[xi].pixelIsInCells = FALSE; /* default - signals that the pixel is outside the cell mesh. */

    }else{
      si = (int)gsl_interp_accel_find(acc, faceDistValues, (size_t)(cellChain->nCellsInChain + 1), x);
      raster[xi].cellAlongRayI = si;
      raster[xi].pixelIsInCells = TRUE;
      cellChain->flags[si] = TRUE;
    }
  }

  gsl_interp_accel_free(acc);
}

/*....................................................................*/
void
interpOnGridAlongRay2ndOrder(baryBuffType *baryBuff, double *vertexValues\
  , double *midEdgeValues, struct simplex *cells, const double deltaX\
  , cellChainType *cellChain, const int numXi, rasterType *raster\
  , double *rasterValues){
  /*
We have a ray which passes through a connected set of simplicial cells. We calculate here the interpolated values for an evenly-spaced sequence of locations along the ray. The coordinate, which we label 'x', is distance along the ray from its origin. The face intercept objects supplied in the arguments entryIntcptFirstCell and cellExitIntcpts already record the value of this coordinate at each intercept point. The number of samples is supplied in argument 'numXi' and their spacing by 'deltaX'.

Input arguments:
================
	- baryBuff: must have been initialized via init_baryBuff(). After that, it contains some useful integers giving the number of spatial dimensions, number of cell vertices etc.

	- vertexValues: contains samples at all cell vertices of the function to be interpolated.

	- midEdgeValues: contains samples at the mid-point of all cell edges (i.e. links between vertices) of the function to be interpolated. Note that such mid-edge values are necessary for 2nd-order interpolation.

	- cells: contains geometrical information about all the cells in the set. Mostly it is used to connect a given cell with the values in 'vertexValues' and 'midEdgeValues' which correspond to that cell's vertices and edges.

	- deltaX: the spacing of the samples.

	- cellChain: contains information (chiefly cell IDs and cell-face ray intercepts) about all the cells traversed by the ray.

	- numXi: required number of samples along the ray.

Output arguments:
=================
	- raster: must have been malloc'd to size numXi. Contains the index (in the lists cellChain->cellIds etc) of the cell each sample falls within, plus a flag to indicate whether it falls within one at all.

	- rasterValues: must have been malloc'd to size numXi*baryBuff->numElements. Contains the output interpolated values.
  */

  double *faceDistValues=NULL,*faceInterpValues=NULL,*midInterpValues=NULL;
  int xi,si,ei,i;
  double x,fracDist,shapeFns[3],ys[3];

  faceDistValues = malloc(sizeof(*faceDistValues)*(cellChain->nCellsInChain+1));

  i = 0;
  faceDistValues[i] = cellChain->entryIntcpt.dist;
  for(i=1;i<=cellChain->nCellsInChain;i++)
    faceDistValues[i] = cellChain->exitIntcpts[i-1].dist;

  midInterpValues  = malloc(sizeof( *midInterpValues)* cellChain->nCellsInChain   *baryBuff->numElements);
  faceInterpValues = malloc(sizeof(*faceInterpValues)*(cellChain->nCellsInChain+1)*baryBuff->numElements);

  _getInterpsAlongRay(baryBuff, vertexValues, midEdgeValues, cells, cellChain, faceInterpValues, midInterpValues);

  _setRasterFlags(baryBuff, deltaX, numXi, faceDistValues, cellChain, raster);

  for(xi=0;xi<numXi;xi++){
    if(!raster[xi].pixelIsInCells){
      for(ei=0;ei<baryBuff->numElements;ei++)
        rasterValues[xi*baryBuff->numElements+ei] = 0.0; /* just to stop it flapping loose. */
  continue;
    }

    x = deltaX*xi;
    si = raster[xi].cellAlongRayI;
    fracDist = (x - faceDistValues[si])/(faceDistValues[si+1] - faceDistValues[si]);
    _getParabolicShapeFns(fracDist, shapeFns);

    for(ei=0;ei<baryBuff->numElements;ei++){
      ys[0] = faceInterpValues[ si   *baryBuff->numElements+ei];
      ys[1] =  midInterpValues[ si   *baryBuff->numElements+ei];
      ys[2] = faceInterpValues[(si+1)*baryBuff->numElements+ei];

      rasterValues[xi*baryBuff->numElements+ei] = interpolateParabolic(ys, shapeFns);
    }
  }

  free(faceDistValues);
  free(faceInterpValues);
  free(midInterpValues);
}



