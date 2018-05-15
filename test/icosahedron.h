#ifndef ICOSAHEDRON_H
#define ICOSAHEDRON_H

#include "../src/raythrucells.h"

int	icosahedronVertices(double **vertexCoords, unsigned long *numPoints);
int	icosahedronGrid(struct gridPoint **gridPoints, unsigned long *numPoints);
int	icosahedronCells(double *vertexCoords, struct simplex **dc, unsigned long *numCells);

#endif /* ICOSAHEDRON_H */

