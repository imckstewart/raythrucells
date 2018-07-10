#include "rtc_types.h"

/*....................................................................*/
struct gridPoint
init_gridPoint(void){
  struct gridPoint gp;

  gp.id = 0;
  gp.numNeigh = 0;
  gp.neigh = NULL;

return gp;
}

/*....................................................................*/
void
free_gridPoints(struct gridPoint *gps, const unsigned long numPoints){
  unsigned long iul;

  if(gps!=NULL){
    for(iul=0;iul<numPoints;iul++)
      free(gps[iul].neigh);

    free(gps);
  }
}

/*....................................................................*/
edgeType
init_edge(void){
  edgeType edge;

  edge.id = 0;
  edge.vertices[0] = 0;
  edge.vertices[1] = 0;

return edge;
}

/*....................................................................*/
struct simplex
init_simplex(void){
  struct simplex cell;
  int i,j,k;

  cell.id = 0;
  i = 0;
  cell.vertx[i] = 0;
  cell.neigh[i] = NULL;
  for(i=1;i<N_DIMS+1;i++){
    cell.vertx[i] = 0;
    cell.neigh[i] = NULL;
    cell.centre[i-1] = 0.0;
  }

  k = 0;
  for(i=0;i<N_DIMS;i++){
    for(j=i+1;j<N_DIMS+1;j++)
      cell.edges[k++] = 0;
  }

return cell;
}

/*....................................................................*/
cellChainType
init_cellChain(const int nCellsMallocd){
  cellChainType chain;

  chain.nCellsInChain = 0;
  chain.nCellsMallocd = nCellsMallocd;
  chain.entryIntcpt = init_intersect();
  if(nCellsMallocd>0){
    chain.exitIntcpts = malloc(sizeof(*chain.exitIntcpts)*nCellsMallocd);
    chain.cellIds     = malloc(sizeof(*chain.cellIds    )*nCellsMallocd);
    chain.flags       = malloc(sizeof(*chain.flags      )*nCellsMallocd);
  }else{
    chain.exitIntcpts = NULL;
    chain.cellIds     = NULL;
    chain.flags       = NULL;
  }

return chain;
}

/*....................................................................*/
cellChainType*
realloc_cellChain(cellChainType *chain, const int nCellsMallocd){
  chain->nCellsMallocd = nCellsMallocd;
  chain->exitIntcpts = realloc(chain->exitIntcpts, sizeof(*chain->exitIntcpts)*nCellsMallocd);
  chain->cellIds     = realloc(chain->cellIds,     sizeof(*chain->cellIds    )*nCellsMallocd);
  chain->flags       = realloc(chain->flags,       sizeof(*chain->flags      )*nCellsMallocd);

return chain;
}

/*....................................................................*/
void
free_cellChain(cellChainType *chain){
  free(chain->exitIntcpts);
  free(chain->cellIds);
  free(chain->flags);
}

/*....................................................................*/
intersectType
init_intersect(void){
  int di;
  intersectType intcpt;

  intcpt.fi = -1;
  intcpt.orientation = 0;
  for(di=0;di<N_DIMS;di++)
    intcpt.bary[di] = 0.0;
  intcpt.dist = 0.0;
  intcpt.collPar = 0.0;

return intcpt;
}

/*....................................................................*/
faceType
init_face(void){
  faceType face;
  int i,j;

  for(i=0;i<N_DIMS;i++){
    face.simplexCentre[i] = 0.0;
    for(j=0;j<N_DIMS;j++)
      face.r[i][j] = 0.0;
  }

return face;
}

/*....................................................................*/
facePlusBasisType
init_facePlusBasis(void){
  facePlusBasisType fpb;
  int i,j;

  for(i=0;i<N_DIMS;i++){
    fpb.origin[i] = 0.0;
    j = 0;
    fpb.axes[i][j] = 0.0;
    for(j=1;j<N_DIMS;j++){
      fpb.axes[i][j] = 0.0;
      fpb.r[i][j-1] = 0.0;
    }
  }

return fpb;
}

/*....................................................................*/
void
_calcEdgeVertexIndices(baryBuffType *baryBuff){
  /*
This should not be called until baryBuff->edgeVertexIndices has been malloc'd, preferably in initializeBaryBuf(). In fact it is intended only to be called from within initializeBaryBuf().
  */
  int ei,i0,i1;

  ei = 0;
  for(i0=0;i0<baryBuff->numVertices-1;i0++){
    for(i1=i0+1;i1<baryBuff->numVertices;i1++){
      baryBuff->edgeVertexIndices[ei][0] = i0;
      baryBuff->edgeVertexIndices[ei][1] = i1;
      ei++;
    }
  }
}

/*....................................................................*/
void
init_baryBuff(const int numDims, const int numElements, baryBuffType *baryBuff){
  /*
This is called only for 2nd-order interpolation.
  */
  baryBuff->numDims = numDims;
  baryBuff->numVertices = numDims + 1;
  baryBuff->numEdges = numDims*(numDims + 1)/2;
  baryBuff->numElements = numElements;
  baryBuff->edgeVertexIndices = malloc(sizeof(*(baryBuff->edgeVertexIndices))*baryBuff->numEdges);
  baryBuff->vertexValues      = malloc(sizeof(*(baryBuff->vertexValues))     *baryBuff->numElements*baryBuff->numVertices);
  baryBuff->edgeValues        = malloc(sizeof(*(baryBuff->edgeValues))       *baryBuff->numElements*baryBuff->numEdges);
  baryBuff->vertxShapeValues  = malloc(sizeof(*(baryBuff->vertxShapeValues)) *baryBuff->numVertices);
  baryBuff->edgeShapeValues   = malloc(sizeof(*(baryBuff->edgeShapeValues))  *baryBuff->numEdges);

  _calcEdgeVertexIndices(baryBuff);
}

/*....................................................................*/
void
free_baryBuff(baryBuffType *baryBuff){
  free(baryBuff->edgeVertexIndices);
  free(baryBuff->vertexValues);
  free(baryBuff->edgeValues);
  free(baryBuff->vertxShapeValues);
  free(baryBuff->edgeShapeValues);
}

