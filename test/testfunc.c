#include "test.h"

_Bool _doTest_tf = FALSE;
int _testFunctionI=0;


/*....................................................................*/
int check_calcEdgeVertexIndices(void){
  int status=0,i;
  const int numDims=3,numElements=1,testEdgeI=2,expectedVertices[]={0,3};
  baryBuffType baryBuff;

  init_baryBuff(numDims, numElements, &baryBuff);

  for(i=0;i<2;i++){
    if(baryBuff.edgeVertexIndices[testEdgeI][i]!=expectedVertices[i]){
      printf("Expected baryBuff.edgeVertexIndices[%d][%d] of %d but got %d\n", testEdgeI, i, expectedVertices[i], baryBuff.edgeVertexIndices[testEdgeI][i]);
      status = 1;
  break;
    }
  }

  free_baryBuff(&baryBuff);

return status;
}


/* Tests of routines in rtc_utils.c:
*/

/*....................................................................*/
int check_gramSchmidt(void){
  int status=0,i;
  const int numDims=3;
  double rawAxes[N_DIMS][N_DIMS],orthoAxes[N_DIMS][N_DIMS];
  double dotProduct;

  rawAxes[0][0] =  3.76;
  rawAxes[0][1] = -1.32;
  rawAxes[0][2] =  2.06;
  rawAxes[1][0] =  1.90;
  rawAxes[1][1] =  1.87;
  rawAxes[1][2] = -0.45;

  gramSchmidt(numDims, numDims-1, rawAxes, orthoAxes);

  /* Check they are orthogonal:
  */
  dotProduct = calcDotProduct(numDims, orthoAxes[0], orthoAxes[1]);

  if(!valuesAreClose(dotProduct, 0.0, 10.0)){
    printf("Axes are not orthogonal - dot product = %e.\n", dotProduct);
return 1;
  }

  /* Check they are of unit length:
  */
  for(i=0;i<numDims-1;i++){
    dotProduct = calcDotProduct(numDims, orthoAxes[i], orthoAxes[i]);

    if(!valuesAreClose(dotProduct, 1.0, 10.0)){
      printf("Axis %d is of length %f, which is not unity.\n", i, sqrt(dotProduct));
return 2;
    }
  }

return status;
}

/*....................................................................*/
int check_calcDotProduct(void){
  int status=0;
  const int numDims=3;
  double dotProduct,expectedDotProduct,vecA[]={0.1,0.5,0.9};

  dotProduct = calcDotProduct(numDims, vecA, vecA);
  expectedDotProduct = 0.01 + 0.25 + 0.81;

  if(!valuesAreClose(dotProduct, expectedDotProduct, 10.0)){
    printf("Expected %f, calcDotProduct() returned %f\n", expectedDotProduct, dotProduct);
return 1;
  }

return status;
}

/*....................................................................*/
int check_calcCellCentres(void){
  int status=0;

  int numDims,numElements;
  struct simplex *cells=NULL;
  unsigned long numCells,numVertices,numEdges;
  double *vertexCoords=NULL,*vertexValues=NULL;
  struct gridPoint *gps=NULL;
  edgeType *edges=NULL;

  int di;
  double expectedCentre[3];

  getExampleCells(&numDims, &numElements, &cells, &numCells, &edges\
    , &numEdges, &vertexCoords, &vertexValues, &gps, &numVertices);

  calcCellCentres(numDims, numCells, vertexCoords, cells);

  for(di=0;di<numDims;di++){
    expectedCentre[di] = (vertexCoords[numDims*0 + di] + vertexCoords[numDims*1 + di]\
                       +  vertexCoords[numDims*2 + di] + vertexCoords[numDims*3 + di])/4.0;

    if(!valuesAreClose(expectedCentre[di], cells[0].centre[di], 10.0)){
      printf("Expected %d coord %f, calcCellCentres() returned %f\n", di, expectedCentre[di], cells[0].centre[di]);
      status = 1;
  break;
    }
  }

  free_gridPoints(gps, numVertices);
  free(edges);
  free(cells);
  free(vertexCoords);
  free(vertexValues);

return status;
}

/*....................................................................*/
int check_getNextEdgeSet(void){
  int status=0,di,i;
  _Bool finished,start=TRUE;
  const int numDims=2,numNeigh=4,expNumSets=6;
  int neighSet[numDims],expectedNeighSets[expNumSets][numDims];

  expectedNeighSets[0][0] = 0;
  expectedNeighSets[0][1] = 1;
  expectedNeighSets[1][0] = 0;
  expectedNeighSets[1][1] = 2;
  expectedNeighSets[2][0] = 0;
  expectedNeighSets[2][1] = 3;
  expectedNeighSets[3][0] = 1;
  expectedNeighSets[3][1] = 2;
  expectedNeighSets[4][0] = 1;
  expectedNeighSets[4][1] = 3;
  expectedNeighSets[5][0] = 2;
  expectedNeighSets[5][1] = 3;

  i = 0;
  while(TRUE){
    finished = _getNextEdgeSet(numDims, numNeigh, &start, neighSet);
  if(finished) break;

    if(i>expNumSets){
      printf("_getNextEdgeSet() did not finish when expected.\n");
return 1;
    }

    for(di=0;di<numDims;di++){
      if(neighSet[di] != expectedNeighSets[i][di]){
        printf("expectedNeighSets[%d][%d] was %d but neighSet[%d] was %d\n", i, di, expectedNeighSets[i][di], di, neighSet[di]);
return 2;
      }
    }

    i++;
  }

  if(i<expNumSets-1){
    printf("_getNextEdgeSet() finished earlier than expected.\n");
return 3;
  }

return status;
}

/*....................................................................*/
int check_gridPointsAreNeighbours(void){
  int status=0;
  const int numDims=3;
  struct gridPoint gpA,gpB,gpC;

  gpA.id = 0;
  gpA.numNeigh = 1;
  gpA.neigh = malloc(sizeof(*(gpA.neigh))*gpA.numNeigh);
  gpA.neigh[0] = &gpB;

  gpB.id = 1;
  gpB.numNeigh = 2;
  gpB.neigh = malloc(sizeof(*(gpB.neigh))*gpB.numNeigh);
  gpB.neigh[0] = &gpA;
  gpB.neigh[1] = &gpC;

  gpC.id = 2;
  gpC.numNeigh = 1;
  gpC.neigh = malloc(sizeof(*(gpC.neigh))*gpC.numNeigh);
  gpC.neigh[0] = &gpB;

  if(!_gridPointsAreNeighbours(numDims, &gpA, &gpB)){
    printf("Grid points A and B should be detected as neighbours.\n");
    status = 1;
  }

  if(status==0 && _gridPointsAreNeighbours(numDims, &gpA, &gpC)){
    printf("Grid points A and C are falsely detected as neighbours.\n");
    status = 2;
  }

  free(gpC.neigh);
  free(gpB.neigh);
  free(gpA.neigh);

return status;
}

/*....................................................................*/
int check_edgesFormACell(void){
  int status=0;

  int numDims,numElements;
  struct simplex *cells=NULL;
  unsigned long numCells,numVertices,numEdges;
  double *vertexCoords=NULL,*vertexValues=NULL;
  struct gridPoint *gps=NULL;
  edgeType *edges=NULL;

  int neighSet[3];

  getExampleCells(&numDims, &numElements, &cells, &numCells, &edges\
    , &numEdges, &vertexCoords, &vertexValues, &gps, &numVertices);

  /* Values refer to indices of gps[3].neigh. */
  neighSet[0] = 0;
  neighSet[1] = 1;
  neighSet[2] = 2;

  if(!_edgesFormACell(numDims, &gps[3], neighSet)){
    printf("Edges should be detected as forming a cell.\n");
    status = 1;
  }

  /* Values refer to indices of gps[3].neigh. */
  neighSet[0] = 0;
  neighSet[1] = 1;
  neighSet[2] = 3;

  if(status==0 && _edgesFormACell(numDims, &gps[3], neighSet)){
    printf("Edges should not be detected as forming a cell.\n");
    status = 2;
  }

  free_gridPoints(gps, numVertices);
  free(edges);
  free(cells);
  free(vertexCoords);
  free(vertexValues);

return status;
}

/*....................................................................*/
int check_cellVerticesMatch(void){
  int status=0;

  int numDims,numElements;
  struct simplex *cells=NULL;
  unsigned long numCells,numVertices,numEdges;
  double *vertexCoords=NULL,*vertexValues=NULL;
  struct gridPoint *gps=NULL;
  edgeType *edges=NULL;

  getExampleCells(&numDims, &numElements, &cells, &numCells, &edges\
    , &numEdges, &vertexCoords, &vertexValues, &gps, &numVertices);

  if(!_cellVerticesMatch(numDims, &cells[0], &cells[0])){
    printf("Cell 0 should match with itself.\n");
    status = 1;
  }

  if(status==0 && _cellVerticesMatch(numDims, &cells[0], &cells[1])){
    printf("Cells 0 and 1 should not match.\n");
    status = 2;
  }

  free_gridPoints(gps, numVertices);
  free(edges);
  free(cells);
  free(vertexCoords);
  free(vertexValues);

return status;
}

/*....................................................................*/
int check_addRawCell(void){
  int status=0;

  int numDims,numElements;
  struct simplex *cells=NULL;
  unsigned long numCells,numVertices,numEdges;
  double *vertexCoords=NULL,*vertexValues=NULL;
  struct gridPoint *gps=NULL;
  edgeType *edges=NULL;

  struct simplex newCell;
  unsigned long maxNumCells=10; /* arbitrary. */

  getExampleCells(&numDims, &numElements, &cells, &numCells, &edges\
    , &numEdges, &vertexCoords, &vertexValues, &gps, &numVertices);

  cells = realloc(cells, sizeof(*cells)*maxNumCells);

  newCell.vertx[0] = 1;
  newCell.vertx[1] = 3;
  newCell.vertx[2] = 4;
  newCell.vertx[3] = 5; /* Non-existent vertex, but this will not matter in the present context. */
  /* The rest of the struct elements are not used in present context. */

  /* Try to add existing cell; should not work. */
  _addRawCell(numDims, &cells[0], &cells, &maxNumCells, &numCells);

  if(numCells>2){
    printf("Cell wrongly added.\n");
    status = 1;
  }

  /* Try to add new cell; should work. */
  _addRawCell(numDims, &newCell, &cells, &maxNumCells, &numCells);

  if(numCells!=3){
    printf("Cell not added.\n");
    status = 2;
  }

  free_gridPoints(gps, numVertices);
  free(edges);
  free(cells);
  free(vertexCoords);
  free(vertexValues);

return status;
}

/*....................................................................*/
int check_getCellsFromGrid(void){
  int status=0;

  int numDims,numElements;
  struct simplex *cells=NULL;
  unsigned long numCells,numVertices,numEdges;
  double *vertexCoords=NULL,*vertexValues=NULL;
  struct gridPoint *gps=NULL;
  edgeType *edges=NULL;

  struct simplex *newCells=NULL;
  unsigned long numNewCells,iul,jul;
  _Bool matchFound;

  getExampleCells(&numDims, &numElements, &cells, &numCells, &edges\
    , &numEdges, &vertexCoords, &vertexValues, &gps, &numVertices);

  getCellsFromGrid(numDims, gps, numVertices, &newCells, &numNewCells);

  if(numNewCells!=numCells){
    printf("Number %lu of new cells should be the same as the number %lu of old.\n", numNewCells, numCells);
    status = 1;
  }

  if(status==0){
    for(iul=0;iul<numCells;iul++){
      matchFound = FALSE; /* default */
      for(jul=0;jul<numNewCells;jul++){
        if(_cellVerticesMatch(numDims, &cells[iul], &newCells[jul])){
          matchFound = TRUE;
      break;
        }
      }
      if(!matchFound){
        printf("Old cell %lu has no match in newCells.\n", iul);
        status = 2;
    break;
      }
    }
  }

  free_gridPoints(gps, numVertices);
  free(edges);
  free(cells);
  free(newCells);
  free(vertexCoords);
  free(vertexValues);

return status;
}

/*....................................................................*/
int check_getEdges(void){
  int status=0;

  int numDims,numElements;
  struct simplex *cells=NULL;
  unsigned long numCells,numVertices,numEdges;
  double *vertexCoords=NULL,*vertexValues=NULL;
  struct gridPoint *gps=NULL;
  edgeType *edges=NULL;

  unsigned long iul;
  int li,ni;
  _Bool otherVertexFound;

  getExampleCells(&numDims, &numElements, &cells, &numCells, &edges\
    , &numEdges, &vertexCoords, &vertexValues, &gps, &numVertices);

  getEdges(numDims, cells, numCells, &edges, &numEdges);

  for(li=0;li<numEdges;li++){
    iul = edges[li].vertices[0];
    otherVertexFound = FALSE; /* default */
    for(ni=0;ni<gps[iul].numNeigh;ni++){
      if(gps[iul].neigh[ni]->id==edges[li].vertices[1]){
        otherVertexFound = TRUE;
    break;
      }
    }

    if(!otherVertexFound){
      printf("For edge %d, could not find a neighbour of point %lu with the expected ID %lu\n", li, iul, edges[li].vertices[1]);
  break;
    }
  }

  free(edges);

  free_gridPoints(gps, numVertices);
  free(edges);
  free(cells);
  free(vertexCoords);
  free(vertexValues);

return status;
}

/*....................................................................*/
int check_calcBaryCoords(void){
  int status=0,i,j;
  const int numDims=3;
  double vertices[numDims+1][numDims],*x=NULL,*bary=NULL;
  double tMat[numDims][numDims],bVec[numDims][1],result[numDims][1],reconstructedX[numDims];

  x    = malloc(sizeof(*x   )* numDims);
  bary = malloc(sizeof(*bary)*(numDims+1));

  vertices[0][0] = 1.0;
  vertices[0][1] = 1.0;
  vertices[0][2] = 0.0;
  vertices[1][0] = 0.0;
  vertices[1][1] = 4.0;
  vertices[1][2] = 0.0;
  vertices[2][0] = 5.0;
  vertices[2][1] = 4.0;
  vertices[2][2] = 0.0;
  vertices[3][0] = 5.0;
  vertices[3][1] = 2.0;
  vertices[3][2] = 1.0;

  x[0] = 2.0;
  x[1] = 2.0;
  x[2] = 0.5;

  calcBaryCoords(numDims, vertices, x, bary);

  /* Back-transform to check them. */
  for(i=0;i<numDims;i++){
    for(j=0;j<numDims;j++)
      tMat[i][j] = vertices[j+1][i] - vertices[0][i];
    bVec[i][0] = bary[i+1];
  }

  matrixMultiply(numDims, 1, numDims, tMat, bVec, result);

  for(i=0;i<numDims;i++){
    reconstructedX[i] = result[i][0] + vertices[0][i];
    if(!valuesAreClose(x[i], reconstructedX[i], 10.0)){
        printf("Reconstructed coord %d is %f, expected to be %f.\n", i, reconstructedX[i], x[i]);
      status = 1;
  break;
    }
  }

  free(bary);
  free(x);

return status;
}



/* Tests of routines in raythrucells.c:
*/

/*....................................................................*/
int check_extractFace(void){
  int status=0;

  int numDims,numElements;
  struct simplex *cells=NULL;
  unsigned long numCells,numVertices,numEdges;
  double *vertexCoords=NULL,*vertexValues=NULL;
  struct gridPoint *gps=NULL;
  edgeType *edges=NULL;

  int vi,vvi,di;
  faceType face;
  unsigned long gi;
  const unsigned long dci=0;
  const int fi=1;

  getExampleCells(&numDims, &numElements, &cells, &numCells, &edges\
    , &numEdges, &vertexCoords, &vertexValues, &gps, &numVertices);

  face = _extractFace(numDims, vertexCoords, cells, dci, fi);

  vvi = 0;
  for(vi=0;vi<numDims+1;vi++){
    if(vi!=fi){
      gi = cells[dci].vertx[vi];
      for(di=0;di<numDims;di++){
        if(!valuesAreClose(vertexCoords[numDims*gi+di], face.r[vvi][di], 10.0)){
          printf("Returned value %f is not close to the expected %f\n", face.r[vvi][di], vertexCoords[numDims*gi+di]);
          status = 1;
      break;
        }
      }
      vvi++;
    }
  if(status!=0) break;
  }

  free_gridPoints(gps, numVertices);
  free(edges);
  free(cells);
  free(vertexCoords);
  free(vertexValues);

return status;
}

/*....................................................................*/
int check_getNewEntryFaceI(void){
  int status=0;

  int numDims,numElements;
  struct simplex *cells=NULL;
  unsigned long numCells,numVertices,numEdges;
  double *vertexCoords=NULL,*vertexValues=NULL;
  struct gridPoint *gps=NULL;
  edgeType *edges=NULL;

  const unsigned long dci=0;
  int newEntryFaceI;
  const int expectedEntryFaceI=3;

  getExampleCells(&numDims, &numElements, &cells, &numCells, &edges\
    , &numEdges, &vertexCoords, &vertexValues, &gps, &numVertices);

  newEntryFaceI = _getNewEntryFaceI(numDims, dci, cells[1]);

  if(newEntryFaceI != expectedEntryFaceI){
    printf("Returned face ID of old cell %d should equal %d\n", newEntryFaceI, expectedEntryFaceI);
    status = 1;
  }

  free_gridPoints(gps, numVertices);
  free(edges);
  free(cells);
  free(vertexCoords);
  free(vertexValues);

return status;
}

/*....................................................................*/
int check_calcFaceInNMinus1(void){
  /* All the numerical values chosen here are arbitrary. */
  int status=0,di,i;
  facePlusBasisType facePlusBasis;
  const int numDims=3,numVertices=3;
  faceType face;
  double axes[N_DIMS-1][N_DIMS],r[N_DIMS][N_DIMS-1],origin[N_DIMS],tempAxis[N_DIMS];
  double dotProduct,vecLength;
  _Bool verbose=FALSE;

  axes[0][0] =  1.0;
  axes[0][1] = -3.0;
  axes[0][2] =  0.5;
  /* Normalize it. */
  dotProduct = calcDotProduct(numDims, axes[0], axes[0]);
  vecLength = sqrt(dotProduct);
  for(di=0;di<numDims;di++)
    axes[0][di] /= vecLength;

  tempAxis[0] = 3.0;
  tempAxis[1] = 1.0;
  tempAxis[2] = 0.5;

  calcCrossProduct(axes[0], tempAxis, axes[1]);
  /* Normalize it. */
  dotProduct = calcDotProduct(numDims, axes[1], axes[1]);
  vecLength = sqrt(dotProduct);
  for(di=0;di<numDims;di++)
    axes[1][di] /= vecLength;

  if(verbose){
    for(i=0;i<2;i++){
      printf("Axis %d:\n", i);
      for(di=0;di<numDims;di++){
        printf("  %+f\n", axes[i][di]);
      }
    }
  }

  /* Check they are orthogonal:
  */
  dotProduct = calcDotProduct(numDims, axes[0], axes[1]);

  if(!valuesAreClose(dotProduct, 0.0, 10.0)){
    printf("Axes are not orthogonal.\n");
return 1;
  }

  /* Check they are of unit length:
  */
  for(i=0;i<numDims-1;i++){
    dotProduct = calcDotProduct(numDims, axes[i], axes[i]);

    if(!valuesAreClose(dotProduct, 1.0, 10.0)){
      printf("Axis %d is of length %f, which is not unity.\n", i, sqrt(dotProduct));
return 2;
    }
  }

  r[0][0] =  0.0;
  r[0][1] =  0.0;
  r[1][0] =  4.2;
  r[1][1] =  0.0;
  r[2][0] =  0.9;
  r[2][1] =  3.5;

  origin[0] =  0.55;
  origin[1] = -0.81;
  origin[2] =  1.06;

  /* Construct the face vertex coordinates:
  */
  for(i=0;i<numVertices;i++){
    for(di=0;di<numDims;di++){
      face.r[i][di] = origin[di] + r[i][0]*axes[0][di] + r[i][1]*axes[1][di];
    }
  }

  facePlusBasis = _calcFaceInNMinus1(numDims, numVertices, &face);

  if(verbose){
    for(i=0;i<2;i++){
      printf("facePlusBasis axis %d:\n", i);
      for(di=0;di<numDims;di++){
        printf("  %+f\n", facePlusBasis.axes[i][di]);
      }
    }
  }

  for(di=0;di<numDims;di++){
    if(!valuesAreClose(facePlusBasis.origin[di], origin[di], 10.0)){
      printf("Origin coord %d is %f, expected to be %f.\n", di, facePlusBasis.origin[di], origin[di]);
return 3;
    }
  }

  for(i=0;i<2;i++){
    for(di=0;di<numDims;di++){
      if(!valuesAreClose(facePlusBasis.axes[i][di], axes[i][di], 10.0)){
        printf("Axis %d coord %d is %f, expected to be %f.\n", i, di, facePlusBasis.axes[i][di], axes[i][di]);
return 4;
      }
    }
  }

  for(i=0;i<numVertices;i++){
    for(di=0;di<numDims-1;di++){
      if(!valuesAreClose(facePlusBasis.r[i][di], r[i][di], 10.0)){
        printf("Vertex %d coord %d is %f, expected to be %f.\n", i, di, facePlusBasis.r[i][di], r[i][di]);
return 5;
      }
    }
  }

return status;
}

/*....................................................................*/
int check_intersectLineWithFace_B(void){
  int status=0,di,j;
  const int numDims=3;
  double extraVertex[numDims],summ;
  double x[]={2.5,3.5,0.0},dir[]={0.0,0.0,1.0};
  faceType face;
  const double epsilon=1.0e-6;
  intersectType intcpt;
  _Bool verbose=TRUE;
  double expectedDist=0.25,expectedCollPar=0.25,expectedBary[]={0.5, 0.25, 0.25};

  face.r[0][0] = 1.0;
  face.r[0][1] = 1.0;
  face.r[0][2] = 0.0;
  face.r[1][0] = 5.0;
  face.r[1][1] = 4.0;
  face.r[1][2] = 0.0;
  face.r[2][0] = 5.0;
  face.r[2][1] = 2.0;
  face.r[2][2] = 1.0;

  extraVertex[0] = 0.0;
  extraVertex[1] = 4.0;
  extraVertex[2] = 0.0;

  for(di=0;di<numDims;di++){
    summ = extraVertex[di];
    for(j=0;j<numDims;j++)
      summ += face.r[j][di];

    face.simplexCentre[di] = 0.25*summ;
  }

  intcpt = _intersectLineWithFace(numDims, x, dir, &face, epsilon);

  if(verbose){
    printf("Intercept:\n  Orientation %d\n  Distance %f\n  Coll. par. %f\n", intcpt.orientation, intcpt.dist, intcpt.collPar);
    printf("  Bary: [%f  %f  %f]\n", intcpt.bary[0], intcpt.bary[1], intcpt.bary[2]);
  }

  if(intcpt.orientation>0){
    printf("Orientation should be <0 (entering cell)\n");
    status = 1;
  }
  if(status==0 && !valuesAreClose(intcpt.dist, expectedDist, 10.0)){
    printf("Intercept distance is %f, expected to be %f.\n", intcpt.dist, expectedDist);
    status = 2;
  }
  if(status==0 && !valuesAreClose(intcpt.collPar, expectedCollPar, 10.0)){
    printf("Intercept collPar is %f, expected to be %f.\n", intcpt.collPar, expectedCollPar);
    status = 3;
  }
  if(status==0){
    for(di=0;di<numDims;di++){
      if(status==0 && !valuesAreClose(intcpt.collPar, expectedCollPar, 10.0)){
        printf("Intercept bary[%d] is %f, expected to be %f.\n", di, intcpt.bary[di], expectedBary[di]);
        status = 4;
    break;
      }
    }
  }

return status;
}

/*....................................................................*/
int check_intersectLineWithFace(void){
  int status=0,di,j;
  const int numDims=3;
  double *x=NULL,*dir=NULL,extraVertex[numDims],summ;
  faceType face;
  const double epsilon=1.0e-6;
  intersectType intcpt;
  _Bool verbose=FALSE;
  double expectedDist=0.25,expectedCollPar=0.25,expectedBary[]={0.5, 0.25, 0.25};

  x   = malloc(sizeof(*x  )*numDims);
  dir = malloc(sizeof(*dir)*numDims);

  x[0] = 1.75;
  x[1] = 2.0;
  x[2] = 0.0;
  dir[0] = 0.0;
  dir[1] = 0.0;
  dir[2] = 1.0;

  face.r[0][0] = 1.0;
  face.r[0][1] = 1.0;
  face.r[0][2] = 0.0;
  face.r[1][0] = 0.0;
  face.r[1][1] = 4.0;
  face.r[1][2] = 0.0;
  face.r[2][0] = 5.0;
  face.r[2][1] = 2.0;
  face.r[2][2] = 1.0;

  extraVertex[0] = 5.0;
  extraVertex[1] = 4.0;
  extraVertex[2] = 0.0;

  for(di=0;di<numDims;di++){
    summ = extraVertex[di];
    for(j=0;j<numDims;j++)
      summ += face.r[j][di];

    face.simplexCentre[di] = 0.25*summ;
  }

  intcpt = _intersectLineWithFace(numDims, x, dir, &face, epsilon);

  if(verbose){
    printf("Intercept:\n  Orientation %d\n  Distance %f\n  Coll. par. %f\n", intcpt.orientation, intcpt.dist, intcpt.collPar);
    printf("  Bary: [%f  %f  %f]\n", intcpt.bary[0], intcpt.bary[1], intcpt.bary[2]);
  }

  if(intcpt.orientation<0){
    printf("Orientation should be >0 (emerging from cell)\n");
    status = 1;
  }
  if(status==0 && !valuesAreClose(intcpt.dist, expectedDist, 10.0)){
    printf("Intercept distance is %f, expected to be %f.\n", intcpt.dist, expectedDist);
    status = 2;
  }
  if(status==0 && !valuesAreClose(intcpt.collPar, expectedCollPar, 10.0)){
    printf("Intercept collPar is %f, expected to be %f.\n", intcpt.collPar, expectedCollPar);
    status = 3;
  }
  if(status==0){
    for(di=0;di<numDims;di++){
      if(status==0 && !valuesAreClose(intcpt.collPar, expectedCollPar, 10.0)){
        printf("Intercept bary[%d] is %f, expected to be %f.\n", di, intcpt.bary[di], expectedBary[di]);
        status = 4;
    break;
      }
    }
  }


  free(dir);
  free(x);

return status;
}

/*....................................................................*/
int check_followGoodChain(void){
  int status=0;

  int numDims,numElements;
  struct simplex *cells=NULL;
  unsigned long numCells,numVertices,numEdges;
  double *vertexCoords=NULL,*vertexValues=NULL;
  struct gridPoint *gps=NULL;
  edgeType *edges=NULL;

  double x[]={2.5,3.5,0.0},dir[]={0.0,0.0,1.0};
  _Bool *cellVisited=NULL;
  unsigned long dci,iul,expectedIds[]={0,1};
  int expectedExitFis[]={0,1},i;
  const double epsilon=1.0e-6;
  cellChainType cellChain;

  const int numFaces=4;
  _Bool chainEndedOk;
  int numMarginalExits,cellEntryFaceI;
  intersectType marginalExitIntcpts[numFaces];

  getExampleCells(&numDims, &numElements, &cells, &numCells, &edges\
    , &numEdges, &vertexCoords, &vertexValues, &gps, &numVertices);

  cellChain = init_cellChain(RTC_BUFFER_SIZE);
  cellVisited = malloc(sizeof(*cellVisited)*numCells);
  for(dci=0;dci<numCells;dci++)
    cellVisited[dci] = FALSE;

  dci = 0; /* Entry cell. */
  cellEntryFaceI = 3;

  chainEndedOk = _followGoodChain(numDims, x, dir, cells, dci, cellEntryFaceI, epsilon\
    , NULL, vertexCoords, &cellVisited, &cellChain, marginalExitIntcpts, &numMarginalExits);

  if(!chainEndedOk){
    printf("Chain did not end correctly.\n");
    status = 1;
  }
  if(status==0 && cellChain.nCellsInChain!=(int)numCells){
    printf("Should be %lu cells in chain but got %d.\n", numCells, cellChain.nCellsInChain);
    status = 2;
  }
  if(status==0 && cellChain.nCellsMallocd!=(int)numCells){
    printf("Cell-chain pointers should have size %lu but got %d.\n", numCells, cellChain.nCellsMallocd);
    status = 3;
  }
  for(i=0;i<cellChain.nCellsInChain;i++){
    iul = cellChain.cellIds[i];
    if(status==0 && iul!=expectedIds[i]){
      printf("Expected cell ID %lu but got %lu.\n", expectedIds[i], iul);
      status = 4;
  break;
    }
    if(status==0 && !cellVisited[iul]){
      printf("Cell %lu not visited.\n", iul);
      status = 5;
  break;
    }
    if(status==0 && cellChain.exitIntcpts[i].fi != expectedExitFis[i]){
      printf("Expected exit-intercept face %d for the %dth cell in the chain but got %d.\n", expectedExitFis[i], i, cellChain.exitIntcpts[i].fi);
      status = 6;
  break;
    }
  }
  if(status==0 && numMarginalExits!=0){
    printf("Shouldn't have any marginal intercepts.\n");
    status = 7;
  }

  free(cellVisited);
  free_cellChain(&cellChain);

  free_gridPoints(gps, numVertices);
  free(edges);
  free(cells);
  free(vertexCoords);
  free(vertexValues);

return status;
}

/*....................................................................*/
int check_buildRayCellChain(void){
  int status=0;

//***

return status;
}

/*....................................................................*/
int check_followRayThroughCells(void){
  int status=0;

  int numDims,numElements;
  struct simplex *cells=NULL;
  unsigned long numCells,numVertices,numEdges;
  double *vertexCoords=NULL,*vertexValues=NULL;
  struct gridPoint *gps=NULL;
  edgeType *edges=NULL;

  double x[]={2.5,3.5,0.0},dir[]={0.0,0.0,1.0};
  const double epsilon=1.0e-6;
  int frtcStatus;
  unsigned long iul,expectedIds[]={0,1};
  int expectedExitFis[]={0,1},i;
  cellChainType cellChain;

  getExampleCells(&numDims, &numElements, &cells, &numCells, &edges\
    , &numEdges, &vertexCoords, &vertexValues, &gps, &numVertices);

  frtcStatus = followRayThroughCells(numDims, x, dir, cells, numCells, epsilon\
    , NULL, vertexCoords, &cellChain);

  if(frtcStatus!=0){
    printf("followRayThroughCells() returned non-zero status %d\n", frtcStatus);
    status = 1;
  }
  if(status==0 && cellChain.nCellsMallocd!=(int)numCells){
    printf("Cell-chain pointers should have size %lu but got %d.\n", numCells, cellChain.nCellsMallocd);
    status = 2;
  }
  for(i=0;i<cellChain.nCellsInChain;i++){
    iul = cellChain.cellIds[i];
    if(status==0 && iul!=expectedIds[i]){
      printf("Expected cell ID %lu but got %lu.\n", expectedIds[i], iul);
      status = 3;
  break;
    }
    if(status==0 && cellChain.exitIntcpts[i].fi != expectedExitFis[i]){
      printf("Expected exit-intercept face %d for the %dth cell in the chain but got %d.\n", expectedExitFis[i], i, cellChain.exitIntcpts[i].fi);
      status = 4;
  break;
    }
  }

  free_cellChain(&cellChain);

  free_gridPoints(gps, numVertices);
  free(edges);
  free(cells);
  free(vertexCoords);
  free(vertexValues);

return status;
}


/* Tests of routines in meshtocube.c:
*/

/*....................................................................*/
int check_interpolateAtFace(void){
  /* All the numerical values chosen here are arbitrary. */
  const int numDims=3,numVertices=10;
  const int numCellVertices=numDims+1;
  intersectType intercept;
  double *vertexValues=NULL,*values=NULL,expectedValue;
  struct simplex cell;
  int i,status=0;

  vertexValues = malloc(sizeof(*vertexValues)*numVertices);
  for(i=0;i<numVertices;i++)
    vertexValues[i] = fmod(sqrt(2.751*i), 1.0);

  for(i=0;i<numCellVertices;i++)
    cell.vertx[i] = 2 + i;

  intercept.bary[0] = 0.145;
  intercept.bary[1] = 0.303;
  intercept.bary[2] = 1.0 - intercept.bary[0] - intercept.bary[1];
  intercept.fi = 2;

  values = malloc(sizeof(*values)*1);

  _interpolateAtFace(numCellVertices, intercept, vertexValues, &cell, values, 1);

  /* Returned value should be equal to:
  */
  expectedValue = intercept.bary[0]*vertexValues[2] + intercept.bary[1]*vertexValues[3] + intercept.bary[2]*vertexValues[5];

  if(!valuesAreClose(expectedValue, values[0], 10.0)){
    printf("Returned value %f is not close to the expected %f\n", values[0], expectedValue);
    status = 1;
  }

  free(values);
  free(vertexValues);

return status;
}

/*....................................................................*/
int check_getFaceInterpsAlongRay(void){
  /* All the numerical values chosen here are arbitrary. */
  int status=0;

  int numDims,numElements;
  struct simplex *cells=NULL;
  unsigned long numCells,numVertices,numEdges;
  double *vertexCoords=NULL,*vertexValues=NULL;
  edgeType *edges=NULL;
  cellChainType cellChain;
  double *rayOrigin=NULL,*rayDir=NULL;

  double *faceInterpValues=NULL,*faceDistValues=NULL;
  double *vvThisFace=NULL,expectedValue;
  int i,j;
  unsigned long dci;
  intersectType intercept;

  getExampleCellsAndRay(&numDims, &numElements, &cells, &numCells, &edges, &numEdges, &vertexCoords\
    , &vertexValues, &numVertices, &cellChain, &rayOrigin, &rayDir);

  faceInterpValues = malloc(sizeof(*faceInterpValues)*(cellChain.nCellsInChain+1)*numElements);
  faceDistValues   = malloc(sizeof(  *faceDistValues)*(cellChain.nCellsInChain+1));

  _getFaceInterpsAlongRay(numDims, vertexValues, cells\
    , &cellChain, faceInterpValues, faceDistValues, numElements);

  /* Just check the first and last values.
  */
  vvThisFace = malloc(sizeof(*vvThisFace)*numDims);

  intercept = cellChain.entryIntcpt;
  dci = 0;

  i = 0;
  for(j=0;j<numDims+1;j++){
  if(j==intercept.fi) continue;
    vvThisFace[i] = vertexValues[cells[dci].vertx[j]];
    i++;
  }
  expectedValue = intercept.bary[0]*vvThisFace[0] + intercept.bary[1]*vvThisFace[1] + intercept.bary[2]*vvThisFace[2];
  if(!valuesAreClose(expectedValue, faceInterpValues[0], 10.0)){
    printf("Returned value faceInterpValues[0] %f is not close to the expected %f\n", faceInterpValues[0], expectedValue);
    status = 1;
  }
  if(!valuesAreClose(intercept.dist, faceDistValues[0], 10.0)){
    printf("Returned value faceDistValues[0] %f is not close to the expected %f\n", faceDistValues[0], intercept.dist);
    status = 2;
  }

  dci = cellChain.nCellsInChain-1;
  intercept = cellChain.exitIntcpts[dci];

  i = 0;
  for(j=0;j<numDims+1;j++){
  if(j==intercept.fi) continue;
    vvThisFace[i] = vertexValues[cells[dci].vertx[j]];
    i++;
  }
  expectedValue = intercept.bary[0]*vvThisFace[0] + intercept.bary[1]*vvThisFace[1] + intercept.bary[2]*vvThisFace[2];
  if(!valuesAreClose(expectedValue, faceInterpValues[cellChain.nCellsInChain], 10.0)){
    printf("Returned value faceInterpValues[-1] %f is not close to the expected %f\n", faceInterpValues[cellChain.nCellsInChain], expectedValue);
    status = 3;
  }
  if(!valuesAreClose(intercept.dist, faceDistValues[cellChain.nCellsInChain], 10.0)){
    printf("Returned value faceDistValues[-1] %f is not close to the expected %f\n", faceDistValues[cellChain.nCellsInChain], intercept.dist);
    status = 4;
  }

  free(vvThisFace);

  free(rayOrigin);
  free(rayDir);
  free(faceDistValues);
  free(faceInterpValues);
  free_cellChain(&cellChain);
  free(vertexValues);
  free(vertexCoords);
  free(cells);
  free(edges);

return status;
}

/*....................................................................*/
int check_interpOnGridAlongRay(void){
  /* All the numerical values chosen here are arbitrary. */
  int status=0;

  int numDims,numElements;
  struct simplex *cells=NULL;
  unsigned long numCells,numVertices,numEdges;
  double *vertexCoords=NULL,*vertexValues=NULL;
  edgeType *edges=NULL;
  cellChainType cellChain;
  double *rayOrigin=NULL,*rayDir=NULL;

  const double deltaX=0.102;
  double *rasterValues=NULL;
  rasterType *raster=NULL;
  const int numSamples=14;
  int startXi,finisXi,si,ei;
  unsigned long dci;
  intersectType startIntcpt,finisIntcpt;
  double startValues[1],finisValues[1],expectedValue;
  double fracDist;

  getExampleCellsAndRay(&numDims, &numElements, &cells, &numCells, &edges, &numEdges, &vertexCoords\
    , &vertexValues, &numVertices, &cellChain, &rayOrigin, &rayDir);

  raster       = malloc(sizeof(*raster)*numSamples);
  rasterValues = malloc(sizeof(*raster)*numSamples);

  _interpOnGridAlongRay(numDims, vertexValues, cells, deltaX, &cellChain, raster\
    , rasterValues, numSamples, numElements);

  /* Ok first let's check that the flagging is sort of ok. I'll find the expected start and end samples and look at their flags and those of their next-outer samples.
  */
  startXi = (int)ceil(cellChain.entryIntcpt.dist/deltaX);
  if(startXi >= numSamples){
    printf("Raster startXi of %d is >= than the raster size %d\n", startXi, numSamples);
    status = 1;
  }
  if(startXi<0) startXi = 0;

  finisXi = (int)floor(cellChain.exitIntcpts[cellChain.nCellsInChain-1].dist/deltaX);
  if(status==0 && finisXi < 0){
    printf("Raster finisXi of %d is less than 0\n", finisXi);
    status = 2;
  }
  if(finisXi>=numSamples) finisXi = numSamples-1;

  if(status==0 && (!raster[startXi].pixelIsInCells || !raster[finisXi].pixelIsInCells)){
    printf("Raster flags unset for startXi %d or finisXi %d\n", startXi, finisXi);
    status = 3;
  }
  if(status==0 && startXi>0 && raster[startXi-1].pixelIsInCells){
    printf("Raster flags set for startXi-1 %d\n", startXi-1);
    status = 4;
  }
  if(status==0 && finisXi<numSamples-1 && raster[finisXi+1].pixelIsInCells){
    printf("Raster flags set for finisXi+1 %d\n", finisXi+1);
    status = 5;
  }

  /* I'll just calculate the expected value of rasterValues[startXi].
  */
  if(status==0){
    si = raster[startXi].cellAlongRayI;
    dci = cellChain.cellIds[si];

    if(si==0){
      startIntcpt = cellChain.entryIntcpt;
      finisIntcpt = cellChain.exitIntcpts[si];
    }else{
      startIntcpt = cellChain.exitIntcpts[si-1];
      finisIntcpt = cellChain.exitIntcpts[si  ];
    }

    _interpolateAtFace(numDims+1, startIntcpt, vertexValues, &cells[dci]\
      , startValues, numElements);

    _interpolateAtFace(numDims+1, finisIntcpt, vertexValues, &cells[dci]\
      , finisValues, numElements);

    fracDist = (startXi*deltaX - startIntcpt.dist)/(finisIntcpt.dist - startIntcpt.dist);

    for(ei=0;ei<numElements;ei++){
      expectedValue =        fracDist *finisValues[ei]\
                    + (1.0 - fracDist)*startValues[ei];

      if(!valuesAreClose(expectedValue, rasterValues[startXi*numElements+ei], 10.0)){
        printf("Returned value rasterValues[%d] %f is not close to the expected %f\n", startXi*numElements+ei, rasterValues[startXi*numElements+ei], expectedValue);
        status = 6;
    break;
      }
    }
  }

  free(rasterValues);
  free(raster);

  free(rayOrigin);
  free(rayDir);
  free_cellChain(&cellChain);
  free(vertexValues);
  free(vertexCoords);
  free(cells);
  free(edges);

return status;
}

/*....................................................................*/
int check_generateVoxelIndex(void){
  int status=0;

  const int numDims=3;
  int xi,yi,zi,pxi[numDims];
  unsigned long ppi,expectedPpi=19;
  axisType axes[N_DIMS];
  _Bool verbose=FALSE;

  if(verbose) printf("Running check_generateVoxelIndex()\n");

  axes[0].numPixels = 4;
  axes[1].numPixels = 3;
  axes[2].numPixels = 2;

  if(verbose){
    for(xi=0;xi<axes[0].numPixels;xi++){
      pxi[0] = xi;
      for(yi=0;yi<axes[1].numPixels;yi++){
        pxi[1] = yi;
        for(zi=0;zi<axes[2].numPixels;zi++){
          pxi[2] = zi;
          ppi = _generateVoxelIndex(numDims, axes, pxi);
          printf("%d %d %d  %2lu\n", xi, yi, zi, ppi);
        }
      }
    }
  }else{
    pxi[0] = 3;
    pxi[1] = 1;
    pxi[2] = 1;

    ppi = _generateVoxelIndex(numDims, axes, pxi);

    if(ppi!=expectedPpi){
      printf("Expected ppi %lu, got %lu\n", expectedPpi, ppi);
return 1;
    }
  }

//***

return status;
}

/*....................................................................*/
int check_generateNextPixelCombo(void){
  int status=0;

  const int numDims=3;
  int zi,pxi[numDims],di,expectedPxi[]={0,2,1};
  unsigned long ppi,expectedPpi=20;
  axisType axes[N_DIMS];
  _Bool finished,verbose=FALSE,limitDims=FALSE;

  if(verbose) printf("Running check_generateNextPixelCombo()\n");

  axes[0].numPixels = 4;
  axes[1].numPixels = 3;
  axes[2].numPixels = 2;

  for(di=0;di<numDims;di++)
    pxi[di] = 0;

  if(verbose){
    ppi = 0;
    finished = FALSE;
    while(!finished){
      if(limitDims){
        for(zi=0;zi<axes[numDims-1].numPixels;zi++){
          printf("%d %d %d  %2lu\n", pxi[0], pxi[1], zi, ppi+zi);
        }

        finished = _generateNextPixelCombo(numDims-1, axes, pxi, &ppi);
        ppi = (unsigned long)axes[numDims-1].numPixels*ppi;
      }else{
        printf("%d %d %d  %2lu\n", pxi[0], pxi[1], pxi[2], ppi);
        finished = _generateNextPixelCombo(numDims, axes, pxi, &ppi);
      }
    }
  }else{
    pxi[0] = 3;
    pxi[1] = 1;
    pxi[2] = 1;

    finished = _generateNextPixelCombo(numDims, axes, pxi, &ppi);

    if(finished){
      printf("_generateNextPixelCombo() should not have finished!\n");
return 1;
    }
    if(ppi!=expectedPpi){
      printf("Expected ppi %lu, got %lu\n", expectedPpi, ppi);
return 2;
    }
    for(di=0;di<numDims;di++){
      if(pxi[di]!=expectedPxi[di]){
        printf("Expected pxi[%d] %d, got %d\n", di, expectedPxi[di], pxi[di]);
return 3;
      }
    }

    for(di=0;di<numDims;di++)
      pxi[di] = axes[di].numPixels - 1;

    finished = _generateNextPixelCombo(numDims, axes, pxi, &ppi);

    if(!finished){
      printf("_generateNextPixelCombo() should have finished!\n");
return 4;
    }
  }

return status;
}

/*....................................................................*/
int check_cellsToHyperCube(void){
  int status=0;

  const int numDims=3;
  double *vertexCoords=NULL,*hypercube=NULL,*vertexValues=NULL;
  unsigned long numPoints=0,numCells=0,ppi,vvi;
  struct simplex *cells=NULL;
  const int numXPixels=100,numYPixels=100,numZPixels=100;
  const double epsilon=1.0e-6;
  double x,y,z,*ray=NULL;
  int numCubes[]={5,5,5},di;
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
    ppi = _generateVoxelIndex(numDims, axes, pxi);
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

        ppi = _generateVoxelIndex(numDims, axes, pxi);
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

return status;
}


/* Tests of routines in second_order.c:
*/

/*....................................................................*/
int check_evaluate2ndOrderShapeFns(void){
  int status=0,i,i0,i1;
  const int numDims=3,numElements=1;
  baryBuffType baryBuff;
  double barys[numDims+1],expectedValue;

  init_baryBuff(numDims, numElements, &baryBuff);

  barys[0] = 0.28;
  barys[1] = 0.11;
  barys[2] = 0.37;
  barys[3] = 1.0 - barys[0] - barys[1] - barys[2];

  _evaluate2ndOrderShapeFns(barys, &baryBuff);

  for(i=0;i<baryBuff.numVertices;i++){
    expectedValue = barys[i]*(2.0*barys[i] - 1.0);

    if(!valuesAreClose(expectedValue, baryBuff.vertxShapeValues[i], 10.0)){
      printf("Returned shape function for vertex %d gives %f but expected %f\n", i, baryBuff.vertxShapeValues[i], expectedValue);
      status = 1;
  break;
    }
  }
  if(status==0){
    for(i=0;i<baryBuff.numEdges;i++){
      i0 = baryBuff.edgeVertexIndices[i][0];
      i1 = baryBuff.edgeVertexIndices[i][1];

      expectedValue = 4.0*barys[i0]*barys[i1];

      if(!valuesAreClose(expectedValue, baryBuff.edgeShapeValues[i], 10.0)){
        printf("Returned shape function for edge %d gives %f but expected %f\n", i, baryBuff.edgeShapeValues[i], expectedValue);
        status = 2;
  break;
      }
    }
  }

  free_baryBuff(&baryBuff);

return status;
}

/*....................................................................*/
double _calcQuadratic(double *x){
  /*
I will assume that we are working with 3 dimensions.

All values are arbitrary.
  */
  const int numDims=3;
  double offsetX[]={1.1,-0.7,2.7},coeffs[]={1.5,0.4,0.8},resultOffset=1.5;
  double delta,result;
  int di;

  result = resultOffset;
  for(di=0;di<numDims;di++){
    delta = x[di] - offsetX[di];
    result += coeffs[di]*delta*delta;
  }

return result;
}

/*....................................................................*/
int check_interpolate2ndOrderCell(void){
  /*
In this test I am going to define a quadratic in 3 dimensions, define a location, then calculate the value of the quadratic at that location. Then I'm going to compare that value to the one obtained from _interpolate2ndOrderCell().
  */
  int status=0,vi,ei,di,i0,i1;
  const int numDims=3,numElements=1;
  double x[]={2.5,3.0,0.2},expectedValue,vertices[numDims+1][numDims],midEdgeLoc[numDims];
  baryBuffType baryBuff;
  double barys[numDims+1],result[1];

  expectedValue = _calcQuadratic(x);

  /* Define the vertex locations.
  */
  vertices[0][0] = 1.0;
  vertices[0][1] = 1.0;
  vertices[0][2] = 0.0;
  vertices[1][0] = 0.0;
  vertices[1][1] = 4.0;
  vertices[1][2] = 0.0;
  vertices[2][0] = 5.0;
  vertices[2][1] = 4.0;
  vertices[2][2] = 0.0;
  vertices[3][0] = 5.0;
  vertices[3][1] = 2.0;
  vertices[3][2] = 1.0;

  /* Set up the buffer to hold the shape functions and vertex values.
  */
  init_baryBuff(numDims, numElements, &baryBuff);

  /* Calculate the vertex and mid-edge values:
  */
  for(vi=0;vi<baryBuff.numVertices;vi++)
    baryBuff.vertexValues[vi] = _calcQuadratic(vertices[vi]);

  for(ei=0;ei<baryBuff.numEdges;ei++){
    i0 = baryBuff.edgeVertexIndices[ei][0];
    i1 = baryBuff.edgeVertexIndices[ei][1];

    for(di=0;di<numDims;di++)
      midEdgeLoc[di] = 0.5*(vertices[i0][di] + vertices[i1][di]);

    baryBuff.edgeValues[ei] = _calcQuadratic(midEdgeLoc);
  }

  calcBaryCoords(numDims, vertices, x, barys);

  _evaluate2ndOrderShapeFns(barys, &baryBuff);

  _interpolate2ndOrderCell(&baryBuff, barys, result);

  if(!valuesAreClose(expectedValue, result[0], 10.0)){
    printf("_interpolate2ndOrderCell() gives %f but expected %f\n", result[0], expectedValue);
    status = 1;
  }

  free_baryBuff(&baryBuff);

return status;
}

/*....................................................................*/
int check_getParabolicShapeFns(void){
  int status=0,i,j;
  const int numValues=3,numShapeFns=3;
  double xs[]={0.0,0.5,1.0},shapeFns[numShapeFns];

  for(i=0;i<numValues;i++){
    _getParabolicShapeFns(xs[i], shapeFns);

    for(j=0;j<numShapeFns;j++){
      if(i==j){
        if(!valuesAreClose(shapeFns[j], 1.0, 10.0)){
          printf("Shape function %d should equal 1 but equals %f\n", j, shapeFns[j]);
          status = 1;
    break;
        }
      }else{
        if(!valuesAreClose(shapeFns[j], 0.0, 10.0)){
          printf("Shape function %d should equal 0 but equals %f\n", j, shapeFns[j]);
          status = 2;
    break;
        }
      }
    }

    if(status!=0)
  break;
  }

return status;
}

/*....................................................................*/
int check_interpolateParabolic(void){
  int status=0,i;
  const int numShapeFns=3;
  const double xOffset=1.5,yOffset=-0.7,coeff=1.2,xLo=-0.2,xHi=1.9;
  double xFrac,x,delta,ys[numShapeFns],expectedResult,shapeFns[numShapeFns],result;

  for(i=0;i<numShapeFns;i++){
    xFrac = 0.5*i;
    x = xLo*(1.0 - xFrac) + xHi*xFrac;
    delta = x - xOffset;
    ys[i] = yOffset + coeff*delta*delta;
  }

  x = 0.38;
  delta = x - xOffset;
  expectedResult = yOffset + coeff*delta*delta;
  xFrac = (x - xLo)/(xHi - xLo);

  _getParabolicShapeFns(xFrac, shapeFns);
  result = interpolateParabolic(ys, shapeFns);

  if(!valuesAreClose(expectedResult, result, 10.0)){
    printf("Expected %f but got %f\n", expectedResult, result);
    status = 1;
  }

return status;
}

/*....................................................................*/
int check_faceBaryToCellBary(void){
  int status=0;

  int numDims,numElements;
  struct simplex *cells=NULL;
  unsigned long numCells,numVertices,numEdges;
  double *vertexCoords=NULL,*vertexValues=NULL;
  edgeType *edges=NULL;
  cellChainType cellChain;
  double *rayOrigin=NULL,*rayDir=NULL;

  const int cellI=0;
  baryBuffType baryBuff;
  double cellBary[N_DIMS+1];
  int fvi,cvi;

  getExampleCellsAndRay(&numDims, &numElements, &cells, &numCells, &edges, &numEdges, &vertexCoords\
    , &vertexValues, &numVertices, &cellChain, &rayOrigin, &rayDir);

  /* Set up the buffer to hold the shape functions and vertex values.
  */
  init_baryBuff(numDims, numElements, &baryBuff);

  _faceBaryToCellBary(&baryBuff, &cellChain.exitIntcpts[cellI], cellBary);

  fvi = 0;
  for(cvi=0;cvi<numDims+1;cvi++){
    if(cvi==cellChain.exitIntcpts[cellI].fi){
      if(cellBary[cvi]!=0.0){
        printf("Expected zero for cell bary %d but got %f\n", cvi, cellBary[cvi]);
        status = 1;
  break;
      }
    }else{
      if(!valuesAreClose(cellChain.exitIntcpts[cellI].bary[fvi], cellBary[cvi], 1.0)){
        printf("Expected %f for cell bary %d but got %f\n", cellChain.exitIntcpts[cellI].bary[fvi], cvi, cellBary[cvi]);
        status = 2;
  break;
      }

      fvi++;
    }
  }

  free_baryBuff(&baryBuff);

  free(rayOrigin);
  free(rayDir);
  free_cellChain(&cellChain);
  free(vertexValues);
  free(vertexCoords);
  free(cells);
  free(edges);

return status;
}

/*....................................................................*/
int check_fillBaryBuffValues(void){
  int status=0;

  baryBuffType baryBuff;
  struct simplex *cells=NULL;
  unsigned long numCells,numVertices,numEdges;
  edgeType *edges=NULL;
  double *vertexCoords=NULL,*vertexValues=NULL,*midEdgeValues=NULL;
  cellChainType cellChain;
  double *rayOrigin=NULL,*rayDir=NULL;

  const unsigned long dci=0;
  const int vi=1,ei=0,li=2;
  int vvi,bvi,eli,bli;

  getExampleCellsAndRay2ndOrder(&baryBuff, &cells, &numCells, &edges\
    , &midEdgeValues, &numEdges, &vertexCoords, &vertexValues, &numVertices\
    , &cellChain, &rayOrigin, &rayDir);

  _fillBaryBuffValues(vertexValues, midEdgeValues, cells, dci, &baryBuff);

  /* Just test one of each:
  */
  vvi = cells[dci].vertx[vi]*baryBuff.numElements + ei;
  bvi = vi*baryBuff.numElements + ei;
  if(!valuesAreClose(vertexValues[vvi], baryBuff.vertexValues[bvi], 10.0)){
    printf("Expected %f for bary vertex value but got %f\n", vertexValues[vvi], baryBuff.vertexValues[bvi]);
    status = 1;
  }
  eli = cells[dci].edges[li]*baryBuff.numElements + ei;
  bli = li*baryBuff.numElements + ei;
  if(status==0 && !valuesAreClose(midEdgeValues[eli], baryBuff.edgeValues[bli], 10.0)){
    printf("Expected %f for bary edge value but got %f\n", midEdgeValues[eli], baryBuff.edgeValues[bli]);
    status = 2;
  }

  free_baryBuff(&baryBuff);

  free(rayOrigin);
  free(rayDir);
  free_cellChain(&cellChain);
  free(midEdgeValues);
  free(vertexValues);
  free(vertexCoords);
  free(cells);
  free(edges);

return status;
}

/*....................................................................*/
int check_getInterpsAlongRay(void){
  int status=0;

  baryBuffType baryBuff;
  struct simplex *cells=NULL;
  unsigned long numCells,numVertices,numEdges;
  edgeType *edges=NULL;
  double *vertexCoords=NULL,*vertexValues=NULL,*midEdgeValues=NULL;
  cellChainType cellChain;
  double *rayOrigin=NULL,*rayDir=NULL;

  double *faceDistValues=NULL,*faceInterpValues=NULL,*midInterpValues=NULL;
  double *interpValues=NULL,barysA[N_DIMS+1],barysB[N_DIMS+1],barysMid[N_DIMS+1];
  const int cellInChainI=0,ei=0,fi=1,mi=0;
  unsigned long dci;
  int vi,i;

  getExampleCellsAndRay2ndOrder(&baryBuff, &cells, &numCells, &edges\
    , &midEdgeValues, &numEdges, &vertexCoords, &vertexValues, &numVertices\
    , &cellChain, &rayOrigin, &rayDir);

  faceDistValues = malloc(sizeof(*faceDistValues)*(cellChain.nCellsInChain+1));

  i = 0;
  faceDistValues[i] = cellChain.entryIntcpt.dist;
  for(i=1;i<cellChain.nCellsInChain+1;i++)
    faceDistValues[i] = cellChain.exitIntcpts[i-1].dist;

  for(i=0;i<cellChain.nCellsInChain;i++)
    cellChain.flags[i] = TRUE; /* just for the purposes of this test. */

  midInterpValues  = malloc(sizeof( *midInterpValues)* cellChain.nCellsInChain   *baryBuff.numElements);
  faceInterpValues = malloc(sizeof(*faceInterpValues)*(cellChain.nCellsInChain+1)*baryBuff.numElements);

  _getInterpsAlongRay(&baryBuff, vertexValues, midEdgeValues, cells, &cellChain\
    , faceInterpValues, midInterpValues);

  /* Just test one of each:
  */
  interpValues = malloc(sizeof(*interpValues)*baryBuff.numElements);
  dci = cellChain.cellIds[cellInChainI];
  _fillBaryBuffValues(vertexValues, midEdgeValues, cells, dci, &baryBuff);

  _faceBaryToCellBary(&baryBuff, &cellChain.entryIntcpt, barysA);

  _faceBaryToCellBary(&baryBuff, &cellChain.exitIntcpts[cellInChainI], barysB);
  _interpolate2ndOrderCell(&baryBuff, barysB, interpValues);

  if(!valuesAreClose(interpValues[ei], faceInterpValues[fi*baryBuff.numElements+ei], 10.0)){
    printf("Expected %f for face interpolated value but got %f\n", interpValues[ei], faceInterpValues[fi*baryBuff.numElements+ei]);
    status = 1;
  }

  for(vi=0;vi<baryBuff.numVertices;vi++)
    barysMid[vi] = 0.5*(barysA[vi] + barysB[vi]);

  _interpolate2ndOrderCell(&baryBuff, barysMid, interpValues);

  if(status==0 && !valuesAreClose(interpValues[ei], midInterpValues[mi*baryBuff.numElements+ei], 10.0)){
    printf("Expected %f for mid interpolated value but got %f\n", interpValues[ei], midInterpValues[mi*baryBuff.numElements+ei]);
    status = 2;
  }

  free(interpValues);

  free(faceInterpValues);
  free(midInterpValues);
  free(faceDistValues);

  free_baryBuff(&baryBuff);

  free(rayOrigin);
  free(rayDir);
  free_cellChain(&cellChain);
  free(midEdgeValues);
  free(vertexValues);
  free(vertexCoords);
  free(cells);
  free(edges);

return status;
}

/*....................................................................*/
int check_setRasterFlags(void){
  int status=0;

//****

return status;
}

/*....................................................................*/
int check_interpOnGridAlongRay2ndOrder(void){
  int status=0;

  baryBuffType baryBuff;
  struct simplex *cells=NULL;
  unsigned long numCells,numVertices,numEdges;
  edgeType *edges=NULL;
  double *vertexCoords=NULL,*vertexValues=NULL,*midEdgeValues=NULL;
  cellChainType cellChain;
  double *rayOrigin=NULL,*rayDir=NULL;

  const double deltaX=0.102;
  double *rasterValues=NULL,barysMid[N_DIMS+1];
  rasterType *raster=NULL;
  const int numSamples=14,ei=0;
  int startXi,finisXi,si,di;
  intersectType startIntcpt,finisIntcpt;
  double expectedValue,fracDist,ys[3],shapeFns[3],interpValues[1];

  getExampleCellsAndRay2ndOrder(&baryBuff, &cells, &numCells, &edges\
    , &midEdgeValues, &numEdges, &vertexCoords, &vertexValues, &numVertices\
    , &cellChain, &rayOrigin, &rayDir);

  raster       = malloc(sizeof(*raster)*numSamples);
  rasterValues = malloc(sizeof(*raster)*numSamples);

  interpOnGridAlongRay2ndOrder(&baryBuff, vertexValues, midEdgeValues\
    , cells, deltaX, &cellChain, numSamples, raster, rasterValues);

  /* Ok first let's check that the flagging is sort of ok. I'll find the expected start and end samples and look at their flags and those of their next-outer samples.
  */
  startXi = (int)ceil(cellChain.entryIntcpt.dist/deltaX);
  if(startXi >= numSamples){
    printf("Raster startXi of %d is >= than the raster size %d\n", startXi, numSamples);
    status = 1;
  }
  if(startXi<0) startXi = 0;

  finisXi = (int)floor(cellChain.exitIntcpts[cellChain.nCellsInChain-1].dist/deltaX);
  if(status==0 && finisXi < 0){
    printf("Raster finisXi of %d is less than 0\n", finisXi);
    status = 2;
  }
  if(finisXi>=numSamples) finisXi = numSamples-1;

  if(status==0 && (!raster[startXi].pixelIsInCells || !raster[finisXi].pixelIsInCells)){
    printf("Raster flags unset for startXi %d or finisXi %d\n", startXi, finisXi);
    status = 3;
  }
  if(status==0 && startXi>0 && raster[startXi-1].pixelIsInCells){
    printf("Raster flags set for startXi-1 %d\n", startXi-1);
    status = 4;
  }
  if(status==0 && finisXi<numSamples-1 && raster[finisXi+1].pixelIsInCells){
    printf("Raster flags set for finisXi+1 %d\n", finisXi+1);
    status = 5;
  }

  /* I'll just calculate the expected value of rasterValues[startXi].
  */
  if(status==0){
    si = raster[startXi].cellAlongRayI;

    if(si==0){
      startIntcpt = cellChain.entryIntcpt;
      finisIntcpt = cellChain.exitIntcpts[si];
    }else{
      startIntcpt = cellChain.exitIntcpts[si-1];
      finisIntcpt = cellChain.exitIntcpts[si  ];
    }

    for(di=0;di<baryBuff.numDims;di++)
      barysMid[di] = 0.5*(startIntcpt.bary[di] + finisIntcpt.bary[di]);

    _interpolate2ndOrderCell(&baryBuff, startIntcpt.bary, interpValues);
    ys[0] = interpValues[ei];
    _interpolate2ndOrderCell(&baryBuff, barysMid, interpValues);
    ys[1] = interpValues[ei];
    _interpolate2ndOrderCell(&baryBuff, finisIntcpt.bary, interpValues);
    ys[2] = interpValues[ei];

    fracDist = (startXi*deltaX - startIntcpt.dist)/(finisIntcpt.dist - startIntcpt.dist);
    _getParabolicShapeFns(fracDist, shapeFns);
    expectedValue = interpolateParabolic(ys, shapeFns);

    if(!valuesAreClose(expectedValue, rasterValues[startXi*baryBuff.numElements+ei], 10.0)){
      printf("Returned value rasterValues[%d] %f is not close to the expected %f\n"\
        , startXi*baryBuff.numElements+ei, rasterValues[startXi*baryBuff.numElements+ei], expectedValue);
      status = 6;
    }
  }

  free(rasterValues);
  free(raster);

  free_baryBuff(&baryBuff);

  free(rayOrigin);
  free(rayDir);
  free_cellChain(&cellChain);
  free(midEdgeValues);
  free(vertexValues);
  free(vertexCoords);
  free(cells);
  free(edges);

return status;
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

  /* Location to write unity-value pixel:
  */
  pxi[0] = 0;
  pxi[1] = 0;
  pxi[2] = numZPixels-1;

  ppi = _generateVoxelIndex(numDims, axes, pxi);
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




