#include "test.h"

int main(int argc, char *argv[]) {
  int testI=0;

  int numDims=3;
  _Bool start,finished;
  int numNeigh,i;
  int neighSet[numDims];

  struct gridPoint *gp=NULL;
  struct simplex candidateCell,*cells=NULL;
  unsigned long maxNumCells,numCells,numGridPoints,iul;

  _Bool isCell;

  struct simplex cellA,cellB;
  _Bool matchFound;


  if(argc>1) testI=atoi(argv[1]);

  switch(testI){
    case 0:
      makeHypercube();
      break;
    case 1:
      checkCubeIndexing2();
      break;
    case 2:
      checkRayInterp();
      break;
    case 3:
      testFits();
      break;
    case 4:
    case 5:
    case 6:
      project(testI-4);
      break;

    /* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .*/
    case 7:
      start = TRUE;
      finished = FALSE;
      numNeigh = 5;

      do{ /* loop over the combinations of N unique neighbours. */
        finished = getNextEdgeSet(numDims, numNeigh, &start, neighSet);
        if(!finished) printf("[%d, %d, %d]\n", neighSet[0],neighSet[1],neighSet[2]);
      }while(!finished);
      break;

    /* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .*/
    case 8:
      numDims = icosahedronGrid(&gp, &numGridPoints);

      neighSet[0] = 0;
      neighSet[1] = 1;
      neighSet[2] = 2;
      isCell = edgesFormACell(numDims, &gp[0], neighSet);
      printf("Is cell? %d (should be 0)\n", (int)isCell);

      neighSet[0] = 0;
      neighSet[1] = 1;
      neighSet[2] = 4;
      isCell = edgesFormACell(numDims, &gp[0], neighSet);
      printf("Is cell? %d (should be 0)\n", (int)isCell);

      neighSet[0] = 0;
      neighSet[1] = 1;
      neighSet[2] = 5;
      isCell = edgesFormACell(numDims, &gp[0], neighSet);
      printf("Is cell? %d (should be 1)\n", (int)isCell);

      neighSet[0] = 0;
      neighSet[1] = 2;
      neighSet[2] = 5;
      isCell = edgesFormACell(numDims, &gp[0], neighSet);
      printf("Is cell? %d (should be 0)\n", (int)isCell);

      free(gp);
      break;

    /* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .*/
    case 9:
      maxNumCells = RTC_BUFFER_SIZE;
      cells = malloc(sizeof(*cells)*maxNumCells);
      numCells = 0;
      numGridPoints = 6;

      gp = malloc(sizeof(*gp)*numGridPoints);
      for(i=0;i<numGridPoints;i++)
        gp[i].id = (unsigned long)i;

      candidateCell.vertx[0] = 1;
      candidateCell.vertx[1] = 4;
      candidateCell.vertx[2] = 2;
      candidateCell.vertx[3] = 3;
      addRawCell(numDims, &candidateCell, &cells, (int *)&maxNumCells, (int *)&numCells);
      printf("numCells=%lu (expect 1)\n", numCells);

      candidateCell.vertx[0] = 2;
      candidateCell.vertx[1] = 3;
      candidateCell.vertx[2] = 4;
      candidateCell.vertx[3] = 5;
      addRawCell(numDims, &candidateCell, &cells, (int *)&maxNumCells, (int *)&numCells);
      printf("numCells=%lu (expect 2)\n", numCells);

      candidateCell.vertx[0] = 4;
      candidateCell.vertx[1] = 1;
      candidateCell.vertx[2] = 3;
      candidateCell.vertx[3] = 2;
      addRawCell(numDims, &candidateCell, &cells, (int *)&maxNumCells, (int *)&numCells);
      printf("numCells=%lu (expect 2)\n", numCells);

      candidateCell.vertx[0] = 2;
      candidateCell.vertx[1] = 0;
      candidateCell.vertx[2] = 3;
      candidateCell.vertx[3] = 1;
      addRawCell(numDims, &candidateCell, &cells, (int *)&maxNumCells, (int *)&numCells);
      printf("numCells=%lu (expect 3)\n", numCells);

      free(gp);
      free(cells);

      break;

    /* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .*/
    case 10:

      cellA.vertx[0] = 0;
      cellA.vertx[1] = 1;
      cellA.vertx[2] = 2;
      cellA.vertx[3] = 3;
      cellB.vertx[0] = 0;
      cellB.vertx[1] = 2;
      cellB.vertx[2] = 1;
      cellB.vertx[3] = 3;
      matchFound = cellVerticesMatch(numDims, &cellA, &cellB);
      printf("matchFound? %d (expect 1)\n", (int)matchFound);

      cellA.vertx[0] = 0;
      cellA.vertx[1] = 4;
      cellA.vertx[2] = 2;
      cellA.vertx[3] = 3;
      cellB.vertx[0] = 0;
      cellB.vertx[1] = 2;
      cellB.vertx[2] = 1;
      cellB.vertx[3] = 3;
      matchFound = cellVerticesMatch(numDims, &cellA, &cellB);
      printf("matchFound? %d (expect 0)\n", (int)matchFound);

      break;

    /* . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .*/
    case 11:
      numDims = icosahedronGrid(&gp, &numGridPoints);
      getCellsFromGrid(numDims, gp, numGridPoints, &cells, &numCells);

      for(iul=0;iul<numCells;iul++){
        printf("Cell %2lu:\n  vertx=[", iul);
        for(i=0;i<numDims;i++){
          printf("%2lu, ", cells[iul].vertx[i]);
        }
        i = numDims;
        printf("%2lu]\n  neigh=[", cells[iul].vertx[i]);
        for(i=0;i<numDims;i++){
          if(cells[iul].neigh[i]==NULL)
            printf(" -, ");
          else
            printf("%2lu, ", cells[iul].neigh[i]->id);
        }
        i = numDims;
        if(cells[iul].neigh[i]==NULL)
          printf(" -]\n");
        else
          printf("%2lu]\n", cells[iul].neigh[i]->id);
      }

      free(cells);
      free(gp);

      break;

    default:
      printf("Test integer %d not recognized.\n", testI);
return 1;
  }

return 0;
}

