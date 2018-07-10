#include "test.h"

/*....................................................................*/
void print_gridPoint(char *spaces, struct gridPoint gp){
  int i;

  printf("%sID = %lu\n", spaces, gp.id);
  printf("%sNum neigh = %d\n", spaces, gp.numNeigh);
  if(gp.neigh==NULL)
    printf("%sNeigh = NULL\n", spaces);
  else{
    for(i=0;i<gp.numNeigh;i++)
      printf("%s  Neigh[%d]->ID = %lu\n", spaces, i, gp.neigh[i]->id);
  }
}

/*....................................................................*/
void print_simplex(char *spaces, struct simplex cell){
/*
struct simplex {
  unsigned long id;
  unsigned long vertx[N_DIMS+1];
  unsigned long edges[N_DIMS*(N_DIMS+1)/2];
  double centre[N_DIMS];
  struct simplex *neigh[N_DIMS+1];
};
*/
}

/*....................................................................*/
void print_intersect(char *spaces, intersectType intcpt){
/*
typedef struct {
  int fi;
  int orientation;
  double bary[N_DIMS], dist, collPar;
} intersectType;
*/
}

/*....................................................................*/
void
_printNeighSet(const int numDims, int neighSet[numDims]){
  int di;

  printf("[");
  for(di=0;di<numDims-1;di++)
    printf("%2d, ", neighSet[di]);
  printf("%2d]\n", neighSet[numDims-1]);
}

/*....................................................................*/
void print_face(char *spaces, const int numDims, faceType face){
/*
typedef struct {
  double r[N_DIMS][N_DIMS], simplexCentre[N_DIMS];
} faceType;
*/

  int i,j;

  for(i=0;i<numDims;i++){
    printf("%sFace r[%d]: [%f", spaces, i, face.r[i][0]);
    for(j=1;j<numDims;j++){
      printf(", %f", face.r[i][j]);
    }
    printf("]\n");
  }

  printf("%ssimplexCentre: [%f", spaces, face.simplexCentre[0]);
  for(j=1;j<numDims;j++){
    printf(", %f", face.simplexCentre[j]);
  }
  printf("]\n");
}

/*....................................................................*/
void
print_cell(const int numDims, const unsigned long iul, struct simplex cell){
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

