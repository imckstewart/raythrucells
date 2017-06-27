#include "test.h"
#include "../src/raythrucells.h"

void rotate(double rotMat[3][3], double inVec[3], double result[3]){
  result[0] = rotMat[0][0]*inVec[0] + rotMat[0][1]*inVec[1] + rotMat[0][2]*inVec[2];
  result[1] = rotMat[1][0]*inVec[0] + rotMat[1][1]*inVec[1] + rotMat[1][2]*inVec[2];
  result[2] = rotMat[2][0]*inVec[0] + rotMat[2][1]*inVec[1] + rotMat[2][2]*inVec[2];
}

int main () {
  const int testI=1;
  int status=0;
  const int numDims=N_DIMS;
  const double scale=1.0,ditherFrac=0.05,pitchDeg=0.0,rollDeg=0.0,yawDeg=0.0;
  const gsl_rng_type *ranNumGenType = gsl_rng_ranlxs2;
  gsl_rng *ran = gsl_rng_alloc(ranNumGenType);
  struct simplex *cells=NULL;
  unsigned long numCells,numPoints;
  double *vertexCoords=NULL,x[numDims],dir[numDims];
  const double epsilon=1.0e-6;
  intersectType entryIntcpt,*cellExitIntcpts=NULL;
  unsigned long *chainOfCellIds=NULL;
  int lenChainPtrs;
  int *sides=NULL,numSides;

  gsl_rng_set(ran, 1237106);

  icosahedron(ditherFrac, ran\
    , &cells, &numCells, &vertexCoords, &numPoints, &sides, &numSides);

  if(testI==0){
    /* Just pass a single ray. */

    x[0] = 0.02;
    x[1] = 0.05;
    x[2] = -2.0;
    dir[0] = 0.0;
    dir[1] = 0.0;
    dir[2] = 1.0;

    status = followRayThroughCells(numDims, x, dir, cells, numCells, epsilon\
      , NULL, vertexCoords, &entryIntcpt, &chainOfCellIds, &cellExitIntcpts\
      , &lenChainPtrs);

    if(status!=0){
      printf("followRayThroughCells status return %d\n", status);
exit(1);
    }

    printf("Number of cells traversed = %lu\n", lenChainPtrs);
    if(lenChainPtrs>0)
      printf("Length of path through cells = %f\n", cellExitIntcpts[lenChainPtrs-1].dist - entryIntcpt.dist);

    free(chainOfCellIds);
    free(cellExitIntcpts);

  }else if(testI==1){
    /*
Pass a grid of rays and plot an image of the path length through the cell mesh of each ray. Actually what I am going to do is make a perspective plot. I'll define a grid of target points centred on the origin, and a vector to the observer, then rotate both. I'll start with the observer located on the +ve X axis, thus the grid in the YZ plane; a +ve yaw rotation is towards the Y, a +ve pitch rotation is towards the Z.
    */
    const double obsDistance=8.0,yawDeg=34.0,ptcDeg=40.0,gridWidth=5.0;
    const int xiSize=400,yiSize=400;
    int xi,yi,device_id,i,j;
    double obs[numDims],cosYaw,sinYaw,cosPtc,sinPtc,rotMat[numDims][numDims];
    double gridPoint[numDims],pathLength,invMat[numDims][numDims];//,*projVertices
    double result[numDims];
    float tr[6],maxValue,worldXLo,worldXHi,worldYLo,worldYHi;
    float *image=malloc(sizeof(double)*xiSize*yiSize);
    float *projVertxXs=NULL,*projVertxYs=NULL,x0,x1,y0,y1;
    char *device = malloc(sizeof(char)*20);
    unsigned long i_ul;

    tr[0] = -gridWidth*(0.5 + 1.5/(float)(xiSize-1));
    tr[1] = (float)gridWidth/(float)(xiSize-1);
    tr[2] = 0.0;
    tr[3] = -gridWidth*(0.5 + 1.5/(float)(yiSize-1));
    tr[4] = 0.0;
    tr[5] = (float)gridWidth/(float)(yiSize-1);

    worldXLo = -gridWidth*0.5;
    worldXHi =  gridWidth*0.5;
    worldYLo = -gridWidth*0.5;
    worldYHi =  gridWidth*0.5;

    sinYaw = sin(yawDeg*PI/180.0);
    cosYaw = cos(yawDeg*PI/180.0);
    sinPtc = sin(ptcDeg*PI/180.0);
    cosPtc = cos(ptcDeg*PI/180.0);

    rotMat[0][0] = cosYaw*cosPtc;
    rotMat[0][1] = -sinYaw;
    rotMat[0][2] = -cosYaw*sinPtc;
    rotMat[1][0] = sinYaw*cosPtc;
    rotMat[1][1] = cosYaw;
    rotMat[1][2] = -sinYaw*sinPtc;
    rotMat[2][0] = sinPtc;
    rotMat[2][1] = 0.0;
    rotMat[2][2] = cosPtc;

    invMat[0][0] = rotMat[0][0];
    invMat[0][1] = rotMat[1][0];
    invMat[0][2] = rotMat[2][0];
    invMat[1][0] = rotMat[0][1];
    invMat[1][1] = rotMat[1][1];
    invMat[1][2] = rotMat[2][1];
    invMat[2][0] = rotMat[0][2];
    invMat[2][1] = rotMat[1][2];
    invMat[2][2] = rotMat[2][2];

    x[0] = obsDistance;
    x[1] = 0.0;
    x[2] = 0.0;

    rotate(rotMat, x, obs);

    projVertxXs = malloc(sizeof(*projVertxXs)*numPoints);
    projVertxYs = malloc(sizeof(*projVertxYs)*numPoints);
    for(i_ul=0;i_ul<numPoints;i_ul++){
      for(j=0;j<numDims;j++)
        x[j] = vertexCoords[numDims*i_ul+j];

      rotate(invMat, x, result);

      projVertxXs[i_ul] = result[1]*obsDistance/(obsDistance - result[0]);
      projVertxYs[i_ul] = result[2]*obsDistance/(obsDistance - result[0]);
    }

    x[0] = 0.0;
    maxValue = 0.0;
    for(xi=0;xi<xiSize;xi++){
      for(yi=0;yi<yiSize;yi++){
        x[1] = gridWidth*((xi/(float)xiSize)-0.5);
        x[2] = gridWidth*((yi/(float)yiSize)-0.5);
        rotate(rotMat, x, gridPoint);

        dir[0] = gridPoint[0] - obs[0];
        dir[1] = gridPoint[1] - obs[1];
        dir[2] = gridPoint[2] - obs[2];

        status = followRayThroughCells(numDims, obs, dir, cells, numCells, epsilon\
          , NULL, vertexCoords, &entryIntcpt, &chainOfCellIds, &cellExitIntcpts\
          , &lenChainPtrs);

        if(status!=0 && status!=RTC_ERR_NO_ENTRIES){
          printf("followRayThroughCells status return %d\n", status);
exit(1);
        }

        if(status==RTC_ERR_NO_ENTRIES)
          pathLength = 0.0;
        else
          pathLength = (float)(cellExitIntcpts[lenChainPtrs-1].dist - entryIntcpt.dist);

        image[yi*xiSize+xi] = pathLength;
        if(pathLength>maxValue) maxValue = pathLength;

        free(chainOfCellIds);
        free(cellExitIntcpts);
      }
    }

    sprintf(device, "/xs");
    device_id = cpgopen(device);
//    cpgsvp(0.1,0.99,0.1,0.99);
    cpgsvp(0.01,0.99,0.01,0.99);
    cpgwnad(worldXLo,worldXHi,worldYLo,worldYHi);
    cpggray(image, xiSize, yiSize, 1, xiSize, 1, yiSize, maxValue, 0.0, tr);

    for(i=0;i<numSides;i++){
      x0 = projVertxXs[sides[2*i  ]];
      x1 = projVertxXs[sides[2*i+1]];
      y0 = projVertxYs[sides[2*i  ]];
      y1 = projVertxYs[sides[2*i+1]];
      cpgmove(x0, y0);
      cpgdraw(x1, y1);
    }

//    cpgbox("BCNT", 0.0, 0, "BCNT", 0.0, 0);
    cpgbox("BC", 0.0, 0, "BC", 0.0, 0);
    cpgclos();

    free(projVertxXs);
    free(projVertxYs);
    free(image);

  }else{
    printf("Unrecognized value %d of the test flag.\n", testI);
exit(1);
  }

  gsl_rng_free(ran);
  free(vertexCoords);
  free(cells);
  free(sides);

  return 0;
}

