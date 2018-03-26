#ifndef RT_UTILS_H
#define RT_UTILS_H

#include "../src/raythrucells.h"

void	calcCellCentres(const int numDims, const unsigned long numCells, double *vertexCoords, struct simplex *dc);
void	getEdges(const int numDims, struct simplex *cells, const unsigned long numCells, edgeType **edges, unsigned long *numEdges);

#endif /* RT_UTILS_H */

