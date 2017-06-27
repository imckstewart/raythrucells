#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>


#include <cpgplot.h>

#include "../src/raythrucells.h"

#define PI                      3.14159265358979323846	/* pi	*/

void icosahedron(const double ditherFrac, gsl_rng *ran\
  , struct simplex **cells, unsigned long *numCells, double **vertexCoords\
  , unsigned long *numPoints, int **sides, int *numSides);
