
#include "rt_utils.h"

/*....................................................................*/
void
calcCellCentres(const int numDims, const unsigned long numCells, double *vertexCoords, struct simplex *dc){
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


