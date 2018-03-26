#ifndef DISSECTED_CUBE_H
#define DISSECTED_CUBE_H

#include "../src/raythrucells.h"

void	dissectedCubeVertices(const int numCubes[3], double **vertexCoords\
  , unsigned long *numPoints);
void	generateDissectedCubeMesh(const int numCubes[3], double *vertexCoords, struct simplex **cells, unsigned long *numCells);

#endif /* DISSECTED_CUBE_H */

