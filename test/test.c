#include "test.h"

int main(int argc, char *argv[]) {
  int testI=0;

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
    default:
      printf("Test integer %d not recognized.\n", testI);
return 1;
  }

return 0;
}

