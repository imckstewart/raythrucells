#include "test.h"

_Bool _doTest_tf = FALSE;
int _testFunctionI=0;

/*....................................................................*/
void matrixMultiply(double a[3][3], double b[3][3], double result[3][3]){
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

  matrixMultiply(yMat, zMat, yzMat);
  matrixMultiply(xMat, yzMat, allMat);

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

  intersectType entryIntcpt,*cellExitIntcpts=NULL;
  unsigned long *chainOfCellIds=NULL;
  int lenChainPtrs;

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
        , numCells, epsilon, NULL, vertexCoords, &entryIntcpt, &chainOfCellIds\
        , &cellExitIntcpts, &lenChainPtrs);

      if(status){
        printf("followRayThroughCells returned with exit status %d\n", status);
exit(1);
      }

      ppi = flattenImageIndices(sizes, xi, yi, 0);
      image[ppi] = 0.0; /* default */
      if(lenChainPtrs>0){
        startZ = entryIntcpt.dist;
        finisZ = cellExitIntcpts[lenChainPtrs-1].dist;
        image[ppi] = finisZ - startZ;
      }

      free(cellExitIntcpts);
      free(chainOfCellIds);
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
void makeHypercube(void){
  const int numDims=3;
  double *vertexCoords=NULL,*hypercube=NULL,*vertexValues=NULL;
  unsigned long numPoints=0,numCells=0,ppi,vvi;
  struct simplex *cells=NULL;
  const int numXPixels=100,numYPixels=100,numZPixels=100;
  const double epsilon=1.0e-6;
  double x,y,z,*ray=NULL;
  int numCubes[]={5,5,5},di,status=0;
  char *fileName0="junk_interp.fits",*fileName1="junk_direct.fits";
  int sizes[numDims];
  double xLos[numDims],xHis[numDims],imageHalfWidth=0.5,imageXOffset=0.5,imageYOffset=0.5,imageZOffset=0.5;
  axisType axes[N_DIMS];

  int pxi[numDims],xi,yi,zi;
  double frac;

  (void)status; /* Stops the 'variable set but not used' warning. */

  printf("Running makeHypercube()\n");

const int testYi=10,testZi=3;

  dissectedCubeVertices(numCubes, &vertexCoords, &numPoints);

  for(ppi=0;ppi<numPoints;ppi++){
    for(di=0;di<numDims;di++){
      vvi = ppi*numDims + di;
      vertexCoords[vvi] /= (double)numCubes[di];
    }
  }

  generateDissectedCubeMesh(numCubes, vertexCoords, &cells, &numCells);

  vertexValues = malloc(sizeof(*vertexValues)*numPoints);
  for(ppi=0;ppi<numPoints;ppi++){
    x = vertexCoords[ppi*numDims + 0];
    y = vertexCoords[ppi*numDims + 1];
    z = vertexCoords[ppi*numDims + 2];
    if(_testFunctionI==0)
      vertexValues[ppi] = testFunctionLinear(x, y, z);
    else if(_testFunctionI==2)
      vertexValues[ppi] = testFunctionGauss(x, y, z);
    else{
      printf("Unrecognized value %d of module-level variable _testFunctionI\n", _testFunctionI);
exit(1);
    }
  }

  axes[0].numPixels = numXPixels;
  axes[1].numPixels = numYPixels;
  axes[2].numPixels = numZPixels;
  for(di=0;di<numDims;di++){
    axes[di].delta = 1.0/(double)axes[di].numPixels;
    axes[di].origin = 0.5*axes[di].delta;
  }

  status = cellsToHyperCube(numDims, vertexCoords, vertexValues, cells, numCells\
    , epsilon, NULL, axes, 1, NULL, &hypercube);

  /* Load and plot a single ray from the hypercube: */
  ray = malloc(sizeof(*ray)*axes[2].numPixels);

  pxi[1] = testYi;
  pxi[2] = testZi;
  for(xi=0;xi<axes[0].numPixels;xi++){
    pxi[0] = xi;
    ppi = generateVoxelIndex(numDims, axes, pxi);
    ray[xi] = hypercube[ppi];
  }
  doTestPlot(ray, axes[0].numPixels);
  free(ray);

  sizes[0] = numXPixels;
  sizes[1] = numYPixels;
  sizes[2] = numZPixels;
  xLos[0] = imageXOffset - imageHalfWidth;
  xHis[0] = imageXOffset + imageHalfWidth;
  xLos[1] = imageYOffset - imageHalfWidth;
  xHis[1] = imageYOffset + imageHalfWidth;
  xLos[2] = imageZOffset - imageHalfWidth;
  xHis[2] = imageZOffset + imageHalfWidth;

  write3Dfits(fileName0, hypercube, sizes, xLos, xHis);

  for(xi=0;xi<axes[0].numPixels;xi++){
    pxi[0] = xi;
    frac = (xi + 0.5)/(double)(axes[0].numPixels);
    x = xLos[0]*(1.0 - frac) + xHis[0]*frac;
    for(yi=0;yi<axes[1].numPixels;yi++){
      pxi[1] = yi;
      frac = (yi + 0.5)/(double)(axes[1].numPixels);
      y = xLos[1]*(1.0 - frac) + xHis[1]*frac;
      for(zi=0;zi<axes[2].numPixels;zi++){
        pxi[2] = zi;
        frac = (zi + 0.5)/(double)(axes[2].numPixels);
        z = xLos[2]*(1.0 - frac) + xHis[2]*frac;

        ppi = generateVoxelIndex(numDims, axes, pxi);
        if(_testFunctionI==0)
          hypercube[ppi] = testFunctionLinear(x, y, z);
        else if(_testFunctionI==2)
          hypercube[ppi] = testFunctionGauss(x, y, z);
        else{
          printf("Unrecognized value %d of module-level variable _testFunctionI\n", _testFunctionI);
exit(1);
        }
      }
    }
  }

  write3Dfits(fileName1, hypercube, sizes, xLos, xHis);

  free(hypercube);
  free(vertexValues);
  free(cells);
  free(vertexCoords);
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
void checkCubeIndexing(void){
  const int numDims=3;
  int xi,yi,zi,pxi[numDims];
  unsigned long ppi;
  axisType axes[N_DIMS];

  printf("Running checkCubeIndexing()\n");

  axes[0].numPixels = 4;
  axes[1].numPixels = 3;
  axes[2].numPixels = 2;

  for(xi=0;xi<axes[0].numPixels;xi++){
    pxi[0] = xi;
    for(yi=0;yi<axes[1].numPixels;yi++){
      pxi[1] = yi;
      for(zi=0;zi<axes[2].numPixels;zi++){
        pxi[2] = zi;
        ppi = generateVoxelIndex(numDims, axes, pxi);
        printf("%d %d %d  %2lu\n", xi, yi, zi, ppi);
      }
    }
  }
}

/*....................................................................*/
void checkCubeIndexing2(void){
  const int numDims=3;
  int zi,pxi[numDims],di;
  unsigned long ppi;
  axisType axes[N_DIMS];
  _Bool finished;

  printf("Running checkCubeIndexing2()\n");

  axes[0].numPixels = 4;
  axes[1].numPixels = 3;
  axes[2].numPixels = 2;

  for(di=0;di<numDims;di++)
    pxi[di] = 0;

  ppi = 0;
  finished = FALSE;
  while(!finished){
    for(zi=0;zi<axes[numDims-1].numPixels;zi++){
      printf("%d %d %d  %2lu\n", pxi[0], pxi[1], zi, ppi+zi);
    }

    finished = _generateNextPixelCombo(numDims-1, axes, pxi, &ppi);
    ppi = (unsigned long)axes[numDims-1].numPixels*ppi;
  }
}

/*....................................................................*/
void checkRayInterp(void){
  const int numDims=3;
  double *vertexCoords=NULL,*vertexValues=NULL;
  unsigned long numPoints=0,numCells=0,ppi,vvi;
  struct simplex *cells=NULL;
  const int numXPixels=100,numYPixels=100,numZPixels=100;
  const double epsilon=1.0e-6;
  double rayDir[]={0.0,0.0,1.0},x,y,z;
  int numCubes[]={5,5,5},di,rtcStatus=0,lenChainPtrs=0;
  axisType axes[N_DIMS];
  int pxi[numDims];
  intersectType entryIntcptFirstCell;
  unsigned long *chainOfCellIds=NULL;
  intersectType *cellExitIntcpts=NULL;
  double rayOrigin[numDims],*rasterValues=NULL;
  rasterType *raster=NULL;

  const int testXi=46,testYi=51;

  printf("Running checkRayInterp()\n");

  dissectedCubeVertices(numCubes, &vertexCoords, &numPoints);

  for(ppi=0;ppi<numPoints;ppi++){
    for(di=0;di<numDims;di++){
      vvi = ppi*numDims + di;
      vertexCoords[vvi] /= (double)numCubes[di];
    }
  }

  generateDissectedCubeMesh(numCubes, vertexCoords, &cells, &numCells);

  vertexValues = malloc(sizeof(*vertexValues)*numPoints);
  for(ppi=0;ppi<numPoints;ppi++){
    x = vertexCoords[ppi*numDims + 0];
    y = vertexCoords[ppi*numDims + 1];
    z = vertexCoords[ppi*numDims + 2];
    if(_testFunctionI==0)
      vertexValues[ppi] = testFunctionLinear(x, y, z);
    else if(_testFunctionI==2)
      vertexValues[ppi] = testFunctionGauss(x, y, z);
    else{
      printf("Unrecognized value %d of module-level variable _testFunctionI\n", _testFunctionI);
exit(1);
    }
  }

  axes[0].numPixels = numXPixels;
  axes[1].numPixels = numYPixels;
  axes[2].numPixels = numZPixels;
  for(di=0;di<numDims;di++){
    axes[di].delta = 1.0/(double)axes[di].numPixels;
    axes[di].origin = 0.5*axes[di].delta;
  }

  for(di=0;di<numDims-1;di++)
    rayDir[di] = 0.0;
  rayDir[numDims-1] = 1.0;

  rayOrigin[numDims-1] = axes[numDims-1].origin; /* for all the rays. */

  pxi[0] = testXi;
  pxi[1] = testYi;
  pxi[2] = 0;

  for(di=0;di<numDims-1;di++)
    rayOrigin[di] = axes[di].origin + pxi[di]*axes[di].delta;

if(1){
  printf("Ray origin XYZ: %e %e %e\n", rayOrigin[0], rayOrigin[1], rayOrigin[2]);
}
  rtcStatus = followRayThroughCells(numDims, rayOrigin, rayDir, cells\
    , numCells, epsilon, NULL, vertexCoords, &entryIntcptFirstCell, &chainOfCellIds\
    , &cellExitIntcpts, &lenChainPtrs);
  if(rtcStatus!=0){
    free(rasterValues);
    free(raster);
    printf("followRayThroughCells returned status %d\n", rtcStatus);
exit(1);
  }

  raster       = malloc(sizeof(*raster)*axes[numDims-1].numPixels);
  rasterValues = malloc(sizeof(*raster)*axes[numDims-1].numPixels);

  _interpolateAlongRay(numDims, vertexValues, cells\
    , axes[numDims-1].delta, entryIntcptFirstCell, chainOfCellIds\
    , cellExitIntcpts, lenChainPtrs, raster\
    , rasterValues, axes[numDims-1].numPixels, 1);

  doTestPlot(rasterValues, axes[numDims-1].numPixels);

  free(rasterValues);
  free(raster);
  free(chainOfCellIds);
  free(cellExitIntcpts);

  free(vertexValues);
  free(cells);
  free(vertexCoords);
}

/*....................................................................*/
void testFits(void){
  /*
I want to check that write3Dfits() produces a cube with the coordinate axes in the correct order. Here I write a cube which is zero everywhere except for a single unity-valued pixel. 
  */
  const int numDims=3;
  char *fileName="junk_1pixel.fits";
  double *hypercube=NULL;
  int sizes[numDims],di,pxi[numDims];
  double xLos[numDims],xHis[numDims],imageHalfWidth=0.5,imageXOffset=0.5,imageYOffset=0.5,imageZOffset=0.5;;
  unsigned long numPointsInCube,ppi;
  const int numXPixels=100,numYPixels=100,numZPixels=100;
  axisType axes[numDims];

  printf("Running testFits()\n");

  axes[0].numPixels = numXPixels;
  axes[1].numPixels = numYPixels;
  axes[2].numPixels = numZPixels;
  for(di=0;di<numDims;di++){
    axes[di].delta = 1.0/(double)axes[di].numPixels;
    axes[di].origin = 0.5*axes[di].delta;
  }

  numPointsInCube = 1;
  for(di=0;di<numDims;di++)
    numPointsInCube *= (unsigned long)axes[di].numPixels;

  hypercube = malloc(sizeof(*hypercube)*numPointsInCube);

  for(ppi=0;ppi<numPointsInCube;ppi++)
    hypercube[ppi] = 0.0;

  // Location to write unity-value pixel:
  pxi[0] = 0;
  pxi[1] = 0;
  pxi[2] = numZPixels-1;

  ppi = generateVoxelIndex(numDims, axes, pxi);
  hypercube[ppi] = 1.0;

  sizes[0] = numXPixels;
  sizes[1] = numYPixels;
  sizes[2] = numZPixels;
  xLos[0] = imageXOffset - imageHalfWidth;
  xHis[0] = imageXOffset + imageHalfWidth;
  xLos[1] = imageYOffset - imageHalfWidth;
  xHis[1] = imageYOffset + imageHalfWidth;
  xLos[2] = imageZOffset - imageHalfWidth;
  xHis[2] = imageZOffset + imageHalfWidth;

  write3Dfits(fileName, hypercube, sizes, xLos, xHis);

  free(hypercube);
}




