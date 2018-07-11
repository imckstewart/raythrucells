#include "test.h"

/*
This contains (or will eventually contain) unit tests for all the routines in the ../src directory.
*/

int main(int argc, char *argv[]) {
  int status=0;
  int moduleI=0,functionI=0;

  if(argc!=3){
    printf("Useage: ./mytest <name of module> <name of function>\n");
return 1;
  }

  if(!getModuleI(argv[1], &moduleI)){
    printf("Module %s not found.\n", argv[1]);
return 2;
  }

  if(!getFunctionI(argv[2], moduleI, &functionI)){
    printf("Function %s of module %s not found.\n", argv[2], argv[1]);
return 3;
  }

  if(moduleI==MOD_RT_UTILS){
    switch(functionI){
      case FUN_GRAMSCHMIDT:
        status = check_gramSchmidt();
        break;

      case FUN_CALCDOTPRODUCT:
        status = check_calcDotProduct();
        break;

      case FUN_CALCCELLCENTRES:
        status = check_calcCellCentres();
        break;

      case FUN_GETNEXTEDGESET:
        status = check_getNextEdgeSet();
        break;

      case FUN_GRIDPOINTSARENEIGHBOURS:
        status = check_gridPointsAreNeighbours();
        break;

      case FUN_EDGESFORMACELL:
        status = check_edgesFormACell();
        break;

      case FUN_CELLVERTICESMATCH:
        status = check_cellVerticesMatch();
        break;

      case FUN_ADDRAWCELL:
        status = check_addRawCell();
        break;

      case FUN_GETCELLSFROMGRID:
        status = check_getCellsFromGrid();
        break;

      case FUN_GETEDGES:
        status = check_getEdges();
        break;

      case FUN_CALCBARYCOORDS:
        status = check_calcBaryCoords();
        break;

      default:
        status = 4;
    }

  }else if(moduleI==MOD_RAYTHRUCELLS){
    switch(functionI){
      case FUN_EXTRACTFACE:
        status = check_extractFace();
        break;

      case FUN_GETNEWENTRYFACEI:
        status = check_getNewEntryFaceI();
        break;

      case FUN_CALCFACEINNMINUS1:
        status = check_calcFaceInNMinus1();
        break;

      case FUN_INTERSECTLINEWITHFACE:
        status = check_intersectLineWithFace();
        break;

      case FUN_FOLLOWGOODCHAIN:
        status = check_followGoodChain();
        break;

      case FUN_BUILDRAYCELLCHAIN:
        status = check_buildRayCellChain();
        break;

      case FUN_FOLLOWRAYTHROUGHCELLS:
        status = check_followRayThroughCells();
        break;

      default:
        status = 4;
    }

  }else if(moduleI==MOD_MESHTOCUBE){
    switch(functionI){
      case FUN_INTERPOLATEATFACE:
        status = check_interpolateAtFace();
        break;

      case FUN_GETFACEINTERPSALONGRAY:
        status = check_getFaceInterpsAlongRay();
        break;

      case FUN_INTERPONGRIDALONGRAY:
        status = check_interpOnGridAlongRay();
        break;

      case FUN_GENERATEVOXELINDEX:
        status = check_generateVoxelIndex();
        break;

      case FUN_GENERATENEXTPIXELCOMBO:
        status = check_generateNextPixelCombo();
        break;

      case FUN_CELLSTOHYPERCUBE:
        status = check_cellsToHyperCube();
        break;

      default:
        status = 4;
    }

  }else if(moduleI==MOD_SECOND_ORDER){
    switch(functionI){
      case FUN_EVALUATE2NDORDERSHAPEFNS:
        status = check_evaluate2ndOrderShapeFns();
        break;

      case FUN_INTERPOLATE2NDORDERCELL:
        status = check_interpolate2ndOrderCell();
        break;

      case FUN_GETPARABOLICSHAPEFNS:
        status = check_getParabolicShapeFns();
        break;

      case FUN_INTERPOLATEPARABOLIC:
        status = check_interpolateParabolic();
        break;

      case FUN_FACEBARYTOCELLBARY:
        status = check_faceBaryToCellBary();
        break;

      case FUN_FILLBARYBUFFVALUES:
        status = check_fillBaryBuffValues();
        break;

      case FUN_GETINTERPSALONGRAY:
        status = check_getInterpsAlongRay();
        break;

      case FUN_SETRASTERFLAGS:
        status = check_setRasterFlags();
        break;

      case FUN_INTERPONGRIDALONGRAY2NDORDER:
        status = check_interpOnGridAlongRay2ndOrder();
        break;

      default:
        status = 4;
    }
  }

return status;
}

