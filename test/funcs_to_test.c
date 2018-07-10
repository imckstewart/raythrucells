#include "funcs_to_test.h"

_Bool
getModuleI(char *moduleName, int *moduleI){
  _Bool nameFound=FALSE;

  *moduleI = 0; /* means none found. */

  if(strcmp(moduleName, "meshtocube")==0){
    *moduleI = MOD_MESHTOCUBE;
return TRUE;
  }

  if(strcmp(moduleName, "raythrucells")==0){
    *moduleI = MOD_RAYTHRUCELLS;
return TRUE;
  }

  if(strcmp(moduleName, "rt_utils")==0){
    *moduleI = MOD_RT_UTILS;
return TRUE;
  }

  if(strcmp(moduleName, "second_order")==0){
    *moduleI = MOD_SECOND_ORDER;
return TRUE;
  }

return nameFound;
}

_Bool
getFunctionI(char *functionName, int moduleI, int *functionI){
  _Bool nameFound=FALSE;

  *functionI = 0; /* means none found. */

  if(moduleI==MOD_RT_UTILS){
    if(strcmp(functionName, "gramSchmidt")==0){
      *functionI = FUN_GRAMSCHMIDT;
return TRUE;
    }

    if(strcmp(functionName, "calcDotProduct")==0){
      *functionI = FUN_CALCDOTPRODUCT;
return TRUE;
    }

    if(strcmp(functionName, "calcCellCentres")==0){
      *functionI = FUN_CALCCELLCENTRES;
return TRUE;
    }

    if(strcmp(functionName, "_getNextEdgeSet")==0){
      *functionI = FUN_GETNEXTEDGESET;
return TRUE;
    }

    if(strcmp(functionName, "_gridPointsAreNeighbours")==0){
      *functionI = FUN_GRIDPOINTSARENEIGHBOURS;
return TRUE;
    }

    if(strcmp(functionName, "_edgesFormACell")==0){
      *functionI = FUN_EDGESFORMACELL;
return TRUE;
    }

    if(strcmp(functionName, "_cellVerticesMatch")==0){
      *functionI = FUN_CELLVERTICESMATCH;
return TRUE;
    }

    if(strcmp(functionName, "_addRawCell")==0){
      *functionI = FUN_ADDRAWCELL;
return TRUE;
    }

    if(strcmp(functionName, "getCellsFromGrid")==0){
      *functionI = FUN_GETCELLSFROMGRID;
return TRUE;
    }

    if(strcmp(functionName, "getEdges")==0){
      *functionI = FUN_GETEDGES;
return TRUE;
    }

    if(strcmp(functionName, "calcBaryCoords")==0){
      *functionI = FUN_CALCBARYCOORDS;
return TRUE;
    }


  }else if(moduleI==MOD_RAYTHRUCELLS){
    if(strcmp(functionName, "_extractFace")==0){
      *functionI = FUN_EXTRACTFACE;
return TRUE;
    }

    if(strcmp(functionName, "_getNewEntryFaceI")==0){
      *functionI = FUN_GETNEWENTRYFACEI;
return TRUE;
    }

    if(strcmp(functionName, "_calcFaceInNMinus1")==0){
      *functionI = FUN_CALCFACEINNMINUS1;
return TRUE;
    }

    if(strcmp(functionName, "_intersectLineWithFace")==0){
      *functionI = FUN_INTERSECTLINEWITHFACE;
return TRUE;
    }

    if(strcmp(functionName, "_followGoodChain")==0){
      *functionI = FUN_FOLLOWGOODCHAIN;
return TRUE;
    }

    if(strcmp(functionName, "_buildRayCellChain")==0){
      *functionI = FUN_BUILDRAYCELLCHAIN;
return TRUE;
    }

    if(strcmp(functionName, "followRayThroughCells")==0){
      *functionI = FUN_FOLLOWRAYTHROUGHCELLS;
return TRUE;
    }


  }else if(moduleI==MOD_MESHTOCUBE){
    if(strcmp(functionName, "_interpolateAtFace")==0){
      *functionI = FUN_INTERPOLATEATFACE;
return TRUE;
    }
    if(strcmp(functionName, "_getFaceInterpsAlongRay")==0){
      *functionI = FUN_GETFACEINTERPSALONGRAY;
return TRUE;
    }
    if(strcmp(functionName, "_interpOnGridAlongRay")==0){
      *functionI = FUN_INTERPONGRIDALONGRAY;
return TRUE;
    }
    if(strcmp(functionName, "_generateVoxelIndex")==0){
      *functionI = FUN_GENERATEVOXELINDEX;
return TRUE;
    }
    if(strcmp(functionName, "_generateNextPixelCombo")==0){
      *functionI = FUN_GENERATENEXTPIXELCOMBO;
return TRUE;
    }
    if(strcmp(functionName, "cellsToHyperCube")==0){
      *functionI = FUN_CELLSTOHYPERCUBE;
return TRUE;
    }


  }else if(moduleI==MOD_SECOND_ORDER){
    if(strcmp(functionName, "_evaluate2ndOrderShapeFns")==0){
      *functionI = FUN_EVALUATE2NDORDERSHAPEFNS;
return TRUE;
    }

    if(strcmp(functionName, "_interpolate2ndOrderCell")==0){
      *functionI = FUN_INTERPOLATE2NDORDERCELL;
return TRUE;
    }

    if(strcmp(functionName, "_getParabolicShapeFns")==0){
      *functionI = FUN_GETPARABOLICSHAPEFNS;
return TRUE;
    }

    if(strcmp(functionName, "interpolateParabolic")==0){
      *functionI = FUN_INTERPOLATEPARABOLIC;
return TRUE;
    }

    if(strcmp(functionName, "_faceBaryToCellBary")==0){
      *functionI = FUN_FACEBARYTOCELLBARY;
return TRUE;
    }

    if(strcmp(functionName, "_fillBaryBuffValues")==0){
      *functionI = FUN_FILLBARYBUFFVALUES;
return TRUE;
    }

    if(strcmp(functionName, "_getInterpsAlongRay")==0){
      *functionI = FUN_GETINTERPSALONGRAY;
return TRUE;
    }

    if(strcmp(functionName, "_setRasterFlags")==0){
      *functionI = FUN_SETRASTERFLAGS;
return TRUE;
    }

    if(strcmp(functionName, "interpOnGridAlongRay2ndOrder")==0){
      *functionI = FUN_INTERPONGRIDALONGRAY2NDORDER;
return TRUE;
    }
  }

return nameFound;
}

