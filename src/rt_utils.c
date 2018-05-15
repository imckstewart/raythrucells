
#include "raythrucells.h"

/*....................................................................*/
void
calcCellCentres(const int numDims, const unsigned long numCells\
  , double *vertexCoords, struct simplex *dc){

  unsigned long i,vi,vvi;
  int j,k;

  for(i=0;i<numCells;i++){
    for(j=0;j<numDims;j++){
      dc[i].centre[j] = 0.0;
      for(k=0;k<numDims+1;k++){
        vi = dc[i].vertx[k];
        vvi = numDims*vi + j;
        dc[i].centre[j] += vertexCoords[vvi];
      }
      dc[i].centre[j] /= (double)(numDims+1);
    }
  }
}

/*....................................................................*/
_Bool getNextEdgeSet(const int numDims, const int numNeigh, _Bool *start\
  , int neighSet[numDims]){
  /*
Example: suppose numDims=3 and numNeigh=5. If start, successive invocations of the routine will write the following to neighSet:

	[0, 1, 2]
	[0, 1, 3]
	[0, 1, 4]
	[0, 2, 3]
	[0, 2, 4]
	[0, 3, 4]
	[1, 2, 3]
	[1, 2, 4]
	[1, 3, 4]
	[2, 3, 4]

I.e., it will eventually return all possible combinations of N separate integers < numNeigh.
  */

  _Bool finished=FALSE,overflow;
  int di,dj;

  if(*start){
    if(numDims>numNeigh)
return TRUE;

    *start = FALSE;
    for(di=0;di<numDims;di++)
      neighSet[di] = di;

  }else{
    finished = TRUE; /* default */

    for(di=numDims-1;di>=0;di--){
      neighSet[di]++;
      if(neighSet[di]>=numNeigh)
    continue;

      overflow = FALSE; /* default */
      for(dj=di+1;dj<numDims;dj++){
        neighSet[dj] = neighSet[dj-1] + 1;
        if(neighSet[dj]>=numNeigh){
          overflow = TRUE;
      break;
        }
      }
      if(overflow)
    continue;

return FALSE;
    }
  }

return finished;
}

void
_printNeighSet(const int numDims, int neighSet[numDims]){
  int di;

  printf("[");
  for(di=0;di<numDims-1;di++)
    printf("%2d, ", neighSet[di]);
  printf("%2d]\n", neighSet[numDims-1]);
}

/*....................................................................*/
_Bool
gridPointsAreNeighbours(const int numDims, struct gridPoint *gpA, struct gridPoint *gpB){
  int i;

  for(i=0;i<gpA->numNeigh;i++){
    if(gpA->neigh[i]->id==gpB->id)
return TRUE;
  }
return FALSE;
}

/*....................................................................*/
_Bool
edgesFormACell(const int numDims, struct gridPoint *gp, int neighSet[numDims]){
  /*
We have here a list neighSet of numDims edges, all of which contain grid point *gp as one of their vertices. I.e. these are lines which radiate out from *gp. The present function determines whether these edges form part of a single simplex. For this to be true, each pair of grid points defined by neighSet should all be neighbours of each other.
  */
  int di,dj;

  for(di=0;di<numDims-1;di++){
    for(dj=di+1;dj<numDims;dj++){
      if(!gridPointsAreNeighbours(numDims, gp->neigh[neighSet[di]], gp->neigh[neighSet[dj]]))
return FALSE;
    }
  }
return TRUE;
}

/*....................................................................*/
_Bool
cellVerticesMatch(const int numDims, struct simplex *cellA, struct simplex *cellB){
  int di,dj;
  _Bool idFound;

  for(di=0;di<numDims+1;di++){
    idFound = FALSE; /* default */
    for(dj=0;dj<numDims+1;dj++){
      if(cellA->vertx[di] == cellB->vertx[dj]){
        idFound = TRUE;
    break;
      }
    }

    if(!idFound)
return FALSE;
  }

return TRUE;
}

/*....................................................................*/
void
addRawCell(const int numDims, struct simplex *candidateCell\
  , struct simplex **cells, int *maxNumCells, int *numCells){
  /*
The purpose of this is to check that the combination of vertices in *pps[1:numDims+1] is not present in any cell in the list *cells. If the set of vertices is found unique, it is stored as a new cell.
  */

  int i;
  _Bool cellAlreadyFound = FALSE; /* default */

  for(i=0;i<*numCells;i++){

    if(cellVerticesMatch(numDims, &(*cells)[i], candidateCell)){
      cellAlreadyFound = TRUE;
  break;
    }
  }

  if(!cellAlreadyFound){
    /* Add the new cell.
    */
    (*numCells)++;
    if(*numCells > *maxNumCells){
      *maxNumCells += RTC_BUFFER_SIZE;
      *cells = realloc(*cells, sizeof(**cells)*(*maxNumCells));
    }

    for(i=0;i<numDims+1;i++){
      (*cells)[*numCells-1].vertx[i] = candidateCell->vertx[i];
      (*cells)[*numCells-1].neigh[i] = NULL;
    }
  }
}

/*....................................................................*/
void
_printCellInfo(const int numDims, const unsigned long iul, struct simplex cell){
  int i;

  printf("Cell %2lu:\n  vertx=[", iul);
  for(i=0;i<numDims;i++){
    printf("%2lu, ", cell.vertx[i]);
  }
  i = numDims;
  printf("%2lu]\n  neigh=[", cell.vertx[i]);
  for(i=0;i<numDims;i++){
    if(cell.neigh[i]==NULL)
      printf(" -, ");
    else
      printf("%2lu, ", cell.neigh[i]->id);
  }
  i = numDims;
  if(cell.neigh[i]==NULL)
    printf(" -]\n");
  else
    printf("%2lu]\n", cell.neigh[i]->id);
}

/*....................................................................*/
void
getCellsFromGrid(const int numDims, struct gridPoint *gp, const unsigned long numPoints\
  , struct simplex **cells, unsigned long *numCells){

  /*
For the sake of brevity I'll represent numDims in these comments by N.

We want here to construct the set of simplicial (i.e. tetrahedral in 3 spatial dimensions) cells that matches the points and links information stored in gp. By 'construct' we mean (i) first of all, calculate the number of cells; (ii) work out which points lie at the N+1 vertices of each cell; (iii) work out which neighbour cells abut the N+1 faces of each cell. The problem is that we will have to deal with reduncancy. The outer loop is over grid points, but since cells have N+1 vertices, looping over the points means we will count each cell N+1 times. Dealing with this will be slow (of order N^2).
  */

  unsigned long maxNumCells,iul,jul;
  int maxNumRawCellsThisPoint,numRawCellsThisPoint,i,di,dj,numCommonVertices;
  struct simplex *cellsThisPoint=NULL,candidateCell;
  _Bool start,finished,*rawCellMatched=NULL;
  _Bool jMatches[numDims+1],vertexIsCommon;
  int neighSet[numDims],diOfUnsharedVertex,djOfUnsharedVertex;

  *numCells = 0;
  maxNumCells = RTC_BUFFER_SIZE;
  *cells = malloc(sizeof(**cells)*maxNumCells);

  maxNumRawCellsThisPoint = RTC_BUFFER_SIZE;
  cellsThisPoint = malloc(sizeof(*cellsThisPoint)*maxNumRawCellsThisPoint);

  for(iul=0;iul<numPoints;iul++){
    numRawCellsThisPoint = 0;

    /* Each cell with point gp[iul] as a vertex will have N edges which have gp[iul] as one of their points. Here we find all sets of N edges starting from gp[iul], then discard those which don't define a cell.
    */
    start = TRUE;
    while(TRUE){ /* loop over the combinations of N unique neighbours. */
      finished = getNextEdgeSet(numDims, gp[iul].numNeigh, &start, neighSet);
      if(finished)
    break;

      if(edgesFormACell(numDims, &gp[iul], neighSet)){
        candidateCell.vertx[0] = gp[iul].id;
        for(di=0;di<numDims;di++)
          candidateCell.vertx[di+1] = gp[iul].neigh[neighSet[di]]->id;

        addRawCell(numDims, &candidateCell, &cellsThisPoint, &maxNumRawCellsThisPoint, &numRawCellsThisPoint);
      }
    } /* end loop over combinations of N unique neighbours. */

    /* Compare each of the cellsThisPoint with the cell list.
    */
    rawCellMatched = malloc(sizeof(*rawCellMatched)*numRawCellsThisPoint);
    for(i=0;i<numRawCellsThisPoint;i++)
      rawCellMatched[i] = FALSE; /* default */

    for(jul=0;jul<*numCells;jul++){
      for(i=0;i<numRawCellsThisPoint;i++){
        if(rawCellMatched[i])
      continue;

        if(cellVerticesMatch(numDims, &cellsThisPoint[i], &(*cells)[jul])){
          rawCellMatched[i] = TRUE;
      break;
        }
      } /* loop over cells around this grid point */
    } /* loop over main list of cells */

    /* Add unmatched grid cells to main list.
    */
    for(i=0;i<numRawCellsThisPoint;i++){
      if(rawCellMatched[i])
    continue;

      (*numCells)++;
      if(*numCells > maxNumCells){
        maxNumCells += RTC_BUFFER_SIZE;
        *cells = realloc(*cells, sizeof(**cells)*maxNumCells);
      }

      for(di=0;di<numDims+1;di++){
        (*cells)[*numCells-1].vertx[di] = cellsThisPoint[i].vertx[di];
        (*cells)[*numCells-1].neigh[di] = NULL;
      }
    }

    free(rawCellMatched);
  } /* Loop over grid points. */

  /* Fill in the cell neighbours:
  */
  for(iul=0;iul<*numCells-1;iul++){
    for(jul=iul+1;jul<*numCells;jul++){
      /* If cells iul and jul share 3 vertices in common, they are neighbours.
      */
      for(dj=0;dj<numDims+1;dj++)
        jMatches[dj] = FALSE;

      numCommonVertices = 0;
      for(di=0;di<numDims+1;di++){
        vertexIsCommon = FALSE; /* default */
        for(dj=0;dj<numDims+1;dj++){
          if((*cells)[iul].vertx[di]==(*cells)[jul].vertx[dj]){
            jMatches[dj] = TRUE;
            numCommonVertices++;
            vertexIsCommon = TRUE;
        break;
          }
        }
        if(!vertexIsCommon)
          diOfUnsharedVertex = di;
      }

      if(numCommonVertices>=numDims){ /* ==numDims should only be possible, not >numDims. == means the cells share a face. */
        for(dj=0;dj<numDims+1;dj++){
          if(!jMatches[dj]){
            djOfUnsharedVertex = dj;
        break;
          }
        }

        (*cells)[iul].neigh[diOfUnsharedVertex] = &(*cells)[jul];
        (*cells)[jul].neigh[djOfUnsharedVertex] = &(*cells)[iul];
      }
    } /* inner loop over cells */
  } /* outer loop over cells */

  free(cellsThisPoint);
}

/*....................................................................*/
void
getEdges(const int numDims, struct simplex *cells, const unsigned long numCells\
  , edgeType **edges, unsigned long *numEdges){

  unsigned long tallyOfEdges=0,ci,ei;
  const int numVertices=numDims+1;
  int i,j,k,vi0,vi1;

  *numEdges = RTC_BUFFER_SIZE;

  *edges = malloc(sizeof(**edges)*(*numEdges));

  tallyOfEdges = 0;
  for(ci=0;ci<numCells;ci++){
    /* Run through each pair of the cell's vertices. */
    k = 0;
    for(i=0;i<numVertices-1;i++){
      vi0 = cells[ci].vertx[i];
      for(j=i+1;j<numVertices;j++){
        vi1 = cells[ci].vertx[j];

        /* Run through all existing edges and see if these 2 points are stored in any. */
        for(ei=0;ei<tallyOfEdges;ei++){
          if(((*edges)[ei].vertices[0]==vi0 && (*edges)[ei].vertices[1]==vi1)\
          || ((*edges)[ei].vertices[0]==vi1 && (*edges)[ei].vertices[1]==vi0)){
        break;
          }
        }

        cells[ci].edges[k++] = ei; /* true whether we found an edge or ran through all the last loop without finding one. */

        if(ei>=tallyOfEdges){ /* need a new edge */
          (*edges)[tallyOfEdges].id = tallyOfEdges;
          (*edges)[tallyOfEdges].vertices[0] = vi0;
          (*edges)[tallyOfEdges].vertices[1] = vi1;
          tallyOfEdges++;

          if(tallyOfEdges>=(*numEdges)){
            *numEdges += RTC_BUFFER_SIZE;
            *edges = realloc(*edges, sizeof(**edges)*(*numEdges));
          }
        }
      }
    }
  }

  *edges = realloc(*edges, sizeof(**edges)*tallyOfEdges);
  *numEdges = tallyOfEdges;

printf("In raythrucells.getEdges(). There are %d edges.\n", (int)tallyOfEdges);
}


