#ifndef ICOSAHEDRON_H
#define ICOSAHEDRON_H

#include "../src/raythrucells.h"

void	icosahedronVertices(double **vertexCoords, unsigned long *numPoints);
//void	icosahedronCells(double *vertexCoords, edgeType **edges, unsigned long *numEdges, struct simplex **dc, unsigned long *numCells);
void	icosahedronCells(double *vertexCoords, struct simplex **dc, unsigned long *numCells);

#endif /* ICOSAHEDRON_H */

