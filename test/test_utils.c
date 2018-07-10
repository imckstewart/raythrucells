#include "test.h"

/*....................................................................*/
double flabs(double inValue){
  /* Returns floating-point absolute value. */
  if(inValue<0.0)
return -inValue;

return inValue;
}

/*....................................................................*/
_Bool valuesAreClose(double valA, double valB, double epsilonMultiplier){
  double absDiff,aveAbs;

  absDiff = flabs(valA - valB);
  aveAbs = 0.5*(flabs(valA) + flabs(valB));
  if(absDiff <= (aveAbs + 1.0)*DBL_EPSILON*epsilonMultiplier)
return TRUE;

return FALSE;
}

/*....................................................................*/
void calcCrossProduct(double *vecA, double *vecB, double *vecResult){
  /* Assumes numDims==3. Doesn't make sense otherwise. */

  vecResult[0] = vecA[1]*vecB[2] - vecA[2]*vecB[1];
  vecResult[1] = vecA[2]*vecB[0] - vecA[0]*vecB[2];
  vecResult[2] = vecA[0]*vecB[1] - vecA[1]*vecB[0];
}

/*....................................................................*/
void matrixMultiply(const int numRowsOut, const int numColsOut\
  , const int numTermsInSum, double a[numRowsOut][numTermsInSum]\
  , double b[numTermsInSum][numColsOut], double result[numRowsOut][numColsOut]){

  int ri,ci,i;

  for(ri=0;ri<numRowsOut;ri++){
    for(ci=0;ci<numColsOut;ci++){
      result[ri][ci] = 0.0;
      for(i=0;i<numTermsInSum;i++)
        result[ri][ci] += a[ri][i]*b[i][ci];
    }
  }
}

/*....................................................................*/
void matrixMultiply3by3(double a[3][3], double b[3][3], double result[3][3]){
  int ri,ci,i;

  for(ri=0;ri<3;ri++){
    for(ci=0;ci<3;ci++){
      result[ri][ci] = 0.0;
      for(i=0;i<3;i++)
        result[ri][ci] += a[ri][i]*b[i][ci];
    }
  }
}

/*....................................................................*/
void rotate(double aboutXDeg, double aboutYDeg, double aboutZDeg, double *vertexCoords, unsigned long numPoints, int numDims){
  double xMat[3][3],yMat[3][3],zMat[3][3],yzMat[3][3],allMat[3][3];
  double cosX,sinX,cosY,sinY,cosZ,sinZ,result;
  unsigned long ppi;
  int di,ddi;
  double convert=M_PI/180.0;

  cosX = cos(aboutXDeg*convert);
  sinX = sin(aboutXDeg*convert);
  cosY = cos(aboutYDeg*convert);
  sinY = sin(aboutYDeg*convert);
  cosZ = cos(aboutZDeg*convert);
  sinZ = sin(aboutZDeg*convert);

  xMat[0][0] =  1.0;
  xMat[0][1] =  0.0;
  xMat[0][2] =  0.0;
  xMat[1][0] =  0.0;
  xMat[1][1] =  cosX;
  xMat[1][2] = -sinX;
  xMat[2][0] =  0.0;
  xMat[2][1] =  sinX;
  xMat[2][2] =  cosX;

  yMat[0][0] =  cosY;
  yMat[0][1] =  0.0;
  yMat[0][2] =  sinY;
  yMat[1][0] =  0.0;
  yMat[1][1] =  1.0;
  yMat[1][2] =  0.0;
  yMat[2][0] = -sinY;
  yMat[2][1] =  0.0;
  yMat[2][2] =  cosY;

  zMat[0][0] =  cosZ;
  zMat[0][1] = -sinZ;
  zMat[0][2] =  0.0;
  zMat[1][0] =  sinZ;
  zMat[1][1] =  cosZ;
  zMat[1][2] =  0.0;
  zMat[2][0] =  0.0;
  zMat[2][1] =  0.0;
  zMat[2][2] =  1.0;

  matrixMultiply3by3(yMat, zMat, yzMat);
  matrixMultiply3by3(xMat, yzMat, allMat);

  for(ppi=0;ppi<numPoints;ppi++){
    for(di=0;di<numDims;di++){
      result = 0.0;
      for(ddi=0;ddi<numDims;ddi++)
        result += allMat[di][ddi]*vertexCoords[numDims*ppi + ddi];
      vertexCoords[numDims*ppi + di] = result;
    }
  }
}

/*....................................................................*/
void project(const int testI){
  double *vertexCoords=NULL,*image=NULL;
  unsigned long numPoints=0,numCells=0;
  struct simplex *cells=NULL;
  const int numDims=3,numXPixels=290,numYPixels=290;
  const double epsilon=1.0e-6;
  double rayStart[numDims],rayDir[]={0.0,0.0,1.0};
  int numCubes[]={5,5,5},status=0;
  char *fileName="junk.fits";
  int sizes[3],xi,yi,vi,i,dummyNumDims;
  double xLos[3],xHis[3],frac,startZ,finisZ,imageHalfWidth,imageXOffset=0.0,imageYOffset=0.0;
  unsigned long ppi;

  cellChainType cellChain;

  (void)dummyNumDims; /* Stops the 'variable set but not used' warning. */

  printf("Running project() with test %d\n", testI);

  switch(testI){
    case 0:
      /* Just a simple tetrahedron.
      */
      printf("  Projecting a tetrahedron.\n");
      imageHalfWidth=2.1;

      numPoints = 4;
      vertexCoords=malloc(sizeof(*vertexCoords)*numDims*numPoints);

      vi = 0;
      vertexCoords[numDims*vi + 0] = -0.8;
      vertexCoords[numDims*vi + 1] = -0.9;
      vertexCoords[numDims*vi + 2] =  0.0;

      vi = 1;
      vertexCoords[numDims*vi + 0] = -0.4;
      vertexCoords[numDims*vi + 1] =  0.7;
      vertexCoords[numDims*vi + 2] =  0.0;

      vi = 2;
      vertexCoords[numDims*vi + 0] =  0.6;
      vertexCoords[numDims*vi + 1] =  0.4;
      vertexCoords[numDims*vi + 2] =  0.0;

      vi = 3;
      vertexCoords[numDims*vi + 0] = -0.3;
      vertexCoords[numDims*vi + 1] =  0.3;
      vertexCoords[numDims*vi + 2] =  0.5;

      numCells = 1;

      cells = malloc(sizeof(*cells)*numCells);

      i = 0;
      cells[i].id = i;
      cells[i].vertx[0] = 0;
      cells[i].vertx[1] = 1;
      cells[i].vertx[2] = 2;
      cells[i].vertx[3] = 3;
      cells[i].neigh[0] = NULL;
      cells[i].neigh[1] = NULL;
      cells[i].neigh[2] = NULL;
      cells[i].neigh[3] = NULL;

      calcCellCentres(numDims, numCells, vertexCoords, cells);

      break;

    case 1:
      /* A rotated icosahedron.
      */
      printf("  Projecting an icosahedron.\n");
      imageHalfWidth=2.1;

      dummyNumDims = icosahedronVertices(&vertexCoords, &numPoints);
      rotate(10.0, 11.0, 12.0, vertexCoords, numPoints, numDims); /* Just a bit of random rotation. */
      dummyNumDims = icosahedronCells(vertexCoords, &cells, &numCells);

      break;
	  
    case 2:
      /* The space here is cuboidal, divided into sub-cuboids by 'numCubes' divisions in the respective dimensions. Each sub-cuboid is further divided into 6 tetrahedral cells. For pure projection purposes we expect the sub-structure to be invisible - the output image should look like a simple x-ray of the main cuboidal space.
      */ 
      printf("  Projecting a cuboidal mesh of tetrahedra.\n");
      imageHalfWidth=3.7;
      imageXOffset = 2.5;
      imageYOffset = 2.47;

      dissectedCubeVertices(numCubes, &vertexCoords, &numPoints);
      printf("Finished calculating vertices.\n");
      rotate(-5.0, 14.0, 9.0, vertexCoords, numPoints, numDims); /* Just a bit of random rotation. */
      printf("Finished rotating.\n");
      generateDissectedCubeMesh(numCubes, vertexCoords, &cells, &numCells);
      printf("Finished calculating mesh.\n");

      break;
	  
    default:
      printf("Test integer %d not recognized.\n", testI);
  }

  sizes[0] = numXPixels;
  sizes[1] = numYPixels;
  sizes[2] = 1;
  xLos[0] = imageXOffset - imageHalfWidth;
  xHis[0] = imageXOffset + imageHalfWidth;
  xLos[1] = imageYOffset - imageHalfWidth;
  xHis[1] = imageYOffset + imageHalfWidth;
  xLos[2] = -1.0;
  xHis[2] = xLos[2];

  image=malloc(sizeof(*image)*numXPixels*numYPixels);
  rayStart[2] = xLos[2];
  for(xi=0;xi<numXPixels;xi++){
    frac = xi/(double)(numXPixels-1);
    rayStart[0] = xLos[0]*(1.0 - frac) + xHis[0]*frac;
    for(yi=0;yi<numYPixels;yi++){
      frac = yi/(double)(numYPixels-1);
      rayStart[1] = xLos[1]*(1.0 - frac) + xHis[1]*frac;

      status = followRayThroughCells(numDims, rayStart, rayDir, cells\
        , numCells, epsilon, NULL, vertexCoords, &cellChain);

      if(status){
        printf("followRayThroughCells returned with exit status %d\n", status);
exit(1);
      }

      ppi = flattenImageIndices(sizes, xi, yi, 0);
      image[ppi] = 0.0; /* default */
      if(cellChain.nCellsInChain>0){
        startZ = cellChain.entryIntcpt.dist;
        finisZ = cellChain.exitIntcpts[cellChain.nCellsInChain-1].dist;
        image[ppi] = finisZ - startZ;
      }

      free_cellChain(&cellChain);
    }
  }

  write3Dfits(fileName, image, sizes, xLos, xHis);

  free(image);
  free(cells);
  free(vertexCoords);
}

/*....................................................................*/
double testFunctionGauss(const double x, const double y, const double z){
  double result;
  const double xSigma=0.3,ySigma=0.2,zSigma=0.43,xOff=0.5,yOff=0.48,zOff=0.51;
  double xFrac,yFrac,zFrac;

  xFrac = (x - xOff)/xSigma;
  yFrac = (y - yOff)/ySigma;
  zFrac = (z - zOff)/zSigma;

  result = exp(-xFrac*xFrac*0.5 - yFrac*yFrac*0.5 - zFrac*zFrac*0.5);

return result;
}

/*....................................................................*/
double testFunctionLinear(const double x, const double y, const double z){
  double result;
  const double xCoeff=0.1,yCoeff=0.2,zCoeff=0.3,xOff=0.5,yOff=0.48,zOff=0.51;

  result = (x - xOff)*xCoeff + (y - yOff)*yCoeff + (z - zOff)*zCoeff;

return result;
}

/*....................................................................*/
void doTestPlot(double *yValues, const int numValues){
  float *xVals=NULL,*yVals=NULL,minYVal=0.0,maxYVal=0.0,yPad;
  int i,devId;

  (void)devId; /* Stops the 'variable set but not used' warning. */

  xVals = malloc(sizeof(*xVals)*numValues);
  yVals = malloc(sizeof(*yVals)*numValues);
  for(i=0;i<numValues;i++){
    xVals[i] = (i+0.5)/(float)numValues;
    yVals[i] = (float)yValues[i];

    if(i==0 || minYVal>yVals[i])
      minYVal = yVals[i];

    if(i==0 || maxYVal<yVals[i])
      maxYVal = yVals[i];
  }

  yPad = 0.05*(maxYVal - minYVal);

  devId = cpgopen("/xs");

  cpgsvp(0.1, 0.999, 0.1, 0.999);
  cpgswin(-0.05, 1.05, minYVal-yPad, maxYVal+yPad);

  cpgpt(numValues, xVals, yVals, -1);

  cpgbox("BCNT", 0.0, 0, "BCNT", 0.0, 0);

  cpgend();

  free(yVals);
  free(xVals);
}

/*....................................................................*/
void
getExampleCells(int *numDims, int *numElements, struct simplex **cells\
  , unsigned long *numCells, edgeType **edges, unsigned long *numEdges\
  , double **vertexCoords, double **vertexValues, struct gridPoint **gps\
  , unsigned long *numVertices){
  /*
All the numerical values chosen here are arbitrary.

The calling routine should free *cells, *vertexCoords, *vertexValues, *gps after use.
 */

  int vi;

  *numDims     = 3;
  *numElements = 1;
  *numVertices = 5;
  *numCells    = 2;

  *cells = malloc(sizeof(**cells)*(*numCells));
  (*cells)[0].vertx[0] = 0;
  (*cells)[0].vertx[1] = 1;
  (*cells)[0].vertx[2] = 2;
  (*cells)[0].vertx[3] = 3;
  (*cells)[1].vertx[0] = 1;
  (*cells)[1].vertx[1] = 2;
  (*cells)[1].vertx[2] = 3;
  (*cells)[1].vertx[3] = 4;

  (*cells)[0].id = 0;
  (*cells)[1].id = 1;

  (*cells)[0].neigh[0] = &(*cells)[1];
  (*cells)[0].neigh[1] = NULL;
  (*cells)[0].neigh[2] = NULL;
  (*cells)[0].neigh[3] = NULL;
  (*cells)[1].neigh[0] = NULL;
  (*cells)[1].neigh[1] = NULL;
  (*cells)[1].neigh[2] = NULL;
  (*cells)[1].neigh[3] = &(*cells)[0];

  getEdges(*numDims, *cells, *numCells, edges, numEdges); /* This also fills (*cells).edges. */

  *vertexCoords = malloc(sizeof(**vertexCoords)*(*numDims)*(*numVertices));

  vi = 0;
  (*vertexCoords)[(*numDims)*vi + 0] = 1.0;
  (*vertexCoords)[(*numDims)*vi + 1] = 1.0;
  (*vertexCoords)[(*numDims)*vi + 2] = 0.0;

  vi = 1;
  (*vertexCoords)[(*numDims)*vi + 0] = 0.0;
  (*vertexCoords)[(*numDims)*vi + 1] = 4.0;
  (*vertexCoords)[(*numDims)*vi + 2] = 0.0;

  vi = 2;
  (*vertexCoords)[(*numDims)*vi + 0] = 5.0;
  (*vertexCoords)[(*numDims)*vi + 1] = 4.0;
  (*vertexCoords)[(*numDims)*vi + 2] = 0.0;

  vi = 3;
  (*vertexCoords)[(*numDims)*vi + 0] = 5.0;
  (*vertexCoords)[(*numDims)*vi + 1] = 2.0;
  (*vertexCoords)[(*numDims)*vi + 2] = 1.0;

  vi = 4;
  (*vertexCoords)[(*numDims)*vi + 0] = 2.0;
  (*vertexCoords)[(*numDims)*vi + 1] = 5.0;
  (*vertexCoords)[(*numDims)*vi + 2] = 2.0;

  calcCellCentres(*numDims, *numCells, *vertexCoords, *cells);

  *vertexValues = malloc(sizeof(**vertexValues)*(*numVertices));
  for(vi=0;vi<(*numVertices);vi++)
    (*vertexValues)[vi] = fmod(sqrt(2.751*vi), 1.0);

  *gps = malloc(sizeof(**gps)*(*numVertices));

  for(vi=0;vi<(*numVertices);vi++)
    (*gps)[vi].id = vi;

  vi = 0;
  (*gps)[vi].numNeigh = 3;
  (*gps)[vi].neigh = malloc(sizeof(*((*gps)[vi].neigh))*(*gps)[vi].numNeigh);
  (*gps)[vi].neigh[0] = &(*gps)[1];
  (*gps)[vi].neigh[1] = &(*gps)[2];
  (*gps)[vi].neigh[2] = &(*gps)[3];

  vi = 1;
  (*gps)[vi].numNeigh = 4;
  (*gps)[vi].neigh = malloc(sizeof(*((*gps)[vi].neigh))*(*gps)[vi].numNeigh);
  (*gps)[vi].neigh[0] = &(*gps)[0];
  (*gps)[vi].neigh[1] = &(*gps)[2];
  (*gps)[vi].neigh[2] = &(*gps)[3];
  (*gps)[vi].neigh[3] = &(*gps)[4];

  vi = 2;
  (*gps)[vi].numNeigh = 4;
  (*gps)[vi].neigh = malloc(sizeof(*((*gps)[vi].neigh))*(*gps)[vi].numNeigh);
  (*gps)[vi].neigh[0] = &(*gps)[0];
  (*gps)[vi].neigh[1] = &(*gps)[1];
  (*gps)[vi].neigh[2] = &(*gps)[3];
  (*gps)[vi].neigh[3] = &(*gps)[4];

  vi = 3;
  (*gps)[vi].numNeigh = 4;
  (*gps)[vi].neigh = malloc(sizeof(*((*gps)[vi].neigh))*(*gps)[vi].numNeigh);
  (*gps)[vi].neigh[0] = &(*gps)[0];
  (*gps)[vi].neigh[1] = &(*gps)[1];
  (*gps)[vi].neigh[2] = &(*gps)[2];
  (*gps)[vi].neigh[3] = &(*gps)[4];

  vi = 4;
  (*gps)[vi].numNeigh = 3;
  (*gps)[vi].neigh = malloc(sizeof(*((*gps)[vi].neigh))*(*gps)[vi].numNeigh);
  (*gps)[vi].neigh[0] = &(*gps)[1];
  (*gps)[vi].neigh[1] = &(*gps)[2];
  (*gps)[vi].neigh[2] = &(*gps)[3];

}

/*....................................................................*/
void
_addRayToCells(int numDims, struct simplex *cells, unsigned long numCells\
  , const double epsilon, double *vertexCoords, cellChainType *cellChain\
  , double **rayOrigin, double **rayDir){

  int rtcStatus=0;

  *rayOrigin = malloc(sizeof(**rayOrigin)*numDims);
  *rayDir    = malloc(sizeof(**rayDir)   *numDims);

  (*rayOrigin)[0] =  2.5;
  (*rayOrigin)[1] =  3.5;
  (*rayOrigin)[2] = -1.0;

  (*rayDir)[0] = 0.0;
  (*rayDir)[1] = 0.0;
  (*rayDir)[2] = 1.0;

  rtcStatus = followRayThroughCells(numDims, *rayOrigin, *rayDir, cells\
    , numCells, epsilon, NULL, vertexCoords, cellChain);

  if(rtcStatus!=0){
    printf("ERROR! followRayThroughCells() exited with non-zero status %d\n", rtcStatus);
exit(1);
  }
}

/*....................................................................*/
void
getExampleCellsAndRay(int *numDims, int *numElements, struct simplex **cells\
  , unsigned long *numCells, edgeType **edges, unsigned long *numEdges\
  , double **vertexCoords, double **vertexValues, unsigned long *numVertices\
  , cellChainType *cellChain, double **rayOrigin, double **rayDir){
  /*
All the numerical values chosen here are arbitrary.

The calling routine should free *cells, *edges, *vertexCoords, *vertexValues, *cellChain, *rayOrigin, *rayDir after use.
 */

  const double epsilon=1.0e-6;
  struct gridPoint *dummyGps=NULL;

  getExampleCells(numDims, numElements, cells, numCells, edges, numEdges\
    , vertexCoords, vertexValues, &dummyGps, numVertices);

  free_gridPoints(dummyGps, *numVertices);

  _addRayToCells(*numDims, *cells, *numCells, epsilon, *vertexCoords, cellChain\
    , rayOrigin, rayDir);
}

/*....................................................................*/
void
getExampleCells2ndOrder(baryBuffType *baryBuff, struct simplex **cells\
  , unsigned long *numCells, edgeType **edges, double **midEdgeValues\
  , unsigned long *numEdges, double **vertexCoords, double **vertexValues\
  , struct gridPoint **gps, unsigned long *numVertices){

  int numDims,numElements;
  unsigned long iul,vi0,vi1;
  double *tempStoreA=NULL,*tempStoreB=NULL,rand;
  int ei;
  const double ditherAmp=0.15;
  const gsl_rng_type *ranNumGenType = gsl_rng_ranlxs2;

  gsl_rng *ran = gsl_rng_alloc(ranNumGenType);
  gsl_rng_set(ran, 1237106) ;

  getExampleCells(&numDims, &numElements, cells, numCells, edges, numEdges\
    , vertexCoords, vertexValues, gps, numVertices);

  tempStoreA = malloc(sizeof(*tempStoreA)*numElements);
  tempStoreB = malloc(sizeof(*tempStoreA)*numElements);

  init_baryBuff(numDims, numElements, baryBuff);

  /* generate midEdgeValues as averages of their vertex values, plus a bit of dither.
  */
  *midEdgeValues = malloc(sizeof(**midEdgeValues)*(*numEdges)*baryBuff->numElements);

  for(iul=0;iul<(*numEdges);iul++){
//printf("Edge %lu (max %lu)\n", iul, *numEdges-1);
    vi0 = (*edges)[iul].vertices[0]; /* for brevity. */
    vi1 = (*edges)[iul].vertices[1]; /* for brevity. */
    for(ei=0;ei<baryBuff->numElements;ei++){
      (*midEdgeValues)[iul*baryBuff->numElements+ei]\
        = 0.5*((*vertexValues)[vi0*baryBuff->numElements+ei]\
        +      (*vertexValues)[vi1*baryBuff->numElements+ei]);
      rand = gsl_rng_uniform(ran);
      (*midEdgeValues)[iul*baryBuff->numElements+ei] *= (1.0 + ditherAmp*(2.0*rand - 1.0));
    }
  }

  free(tempStoreB);
  free(tempStoreA);

  gsl_rng_free(ran);
}

/*....................................................................*/
void
getExampleCellsAndRay2ndOrder(baryBuffType *baryBuff, struct simplex **cells\
  , unsigned long *numCells, edgeType **edges, double **midEdgeValues\
  , unsigned long *numEdges, double **vertexCoords, double **vertexValues\
  , unsigned long *numVertices, cellChainType *cellChain\
  , double **rayOrigin, double **rayDir){

  const double epsilon=1.0e-6;
  struct gridPoint *dummyGps=NULL;

  getExampleCells2ndOrder(baryBuff, cells, numCells, edges, midEdgeValues\
    , numEdges, vertexCoords, vertexValues, &dummyGps, numVertices);

  free_gridPoints(dummyGps, *numVertices);

  _addRayToCells(baryBuff->numDims, *cells, *numCells, epsilon, *vertexCoords, cellChain\
    , rayOrigin, rayDir);
}



