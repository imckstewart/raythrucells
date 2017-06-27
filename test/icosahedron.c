#include <gsl/gsl_rng.h>
#include "../src/raythrucells.h"

void icosahedron(const double ditherFrac, gsl_rng *randGen\
  , struct simplex **cells, unsigned long *numCells, double **vertexCoords\
  , unsigned long *numPoints, int **sides, int *numSides){
  /*
This constructs vertices and simplicial cells enough to define an icosahedron.
  */
  const int numDims=3;
  int i,j,k,vi;
  const double phi=0.5*(1.0 + sqrt(5.0));
  const double oneOnNumVertices = 1.0/(double)(numDims + 1);
  unsigned long i_ul;

  *numCells = 20;
  *numPoints = 13;

  *vertexCoords = malloc(sizeof(**vertexCoords)*numDims*(*numPoints));

  i = 0;
  (*vertexCoords)[numDims*i+0] =  0.0;
  (*vertexCoords)[numDims*i+1] =  0.0;
  (*vertexCoords)[numDims*i+2] =  0.0;
  i++;
  (*vertexCoords)[numDims*i+0] =  0.0;
  (*vertexCoords)[numDims*i+1] =  1.0;
  (*vertexCoords)[numDims*i+2] =  phi;
  i++;
  (*vertexCoords)[numDims*i+0] =  0.0;
  (*vertexCoords)[numDims*i+1] =  1.0;
  (*vertexCoords)[numDims*i+2] = -phi;
  i++;
  (*vertexCoords)[numDims*i+0] =  0.0;
  (*vertexCoords)[numDims*i+1] = -1.0;
  (*vertexCoords)[numDims*i+2] =  phi;
  i++;
  (*vertexCoords)[numDims*i+0] =  0.0;
  (*vertexCoords)[numDims*i+1] = -1.0;
  (*vertexCoords)[numDims*i+2] = -phi;
  i++;
  (*vertexCoords)[numDims*i+0] =  1.0;
  (*vertexCoords)[numDims*i+1] =  phi;
  (*vertexCoords)[numDims*i+2] =  0.0;
  i++;
  (*vertexCoords)[numDims*i+0] =  1.0;
  (*vertexCoords)[numDims*i+1] = -phi;
  (*vertexCoords)[numDims*i+2] =  0.0;
  i++;
  (*vertexCoords)[numDims*i+0] = -1.0;
  (*vertexCoords)[numDims*i+1] =  phi;
  (*vertexCoords)[numDims*i+2] =  0.0;
  i++;
  (*vertexCoords)[numDims*i+0] = -1.0;
  (*vertexCoords)[numDims*i+1] = -phi;
  (*vertexCoords)[numDims*i+2] =  0.0;
  i++;
  (*vertexCoords)[numDims*i+0] =  phi;
  (*vertexCoords)[numDims*i+1] =  0.0;
  (*vertexCoords)[numDims*i+2] =  1.0;
  i++;
  (*vertexCoords)[numDims*i+0] =  phi;
  (*vertexCoords)[numDims*i+1] =  0.0;
  (*vertexCoords)[numDims*i+2] = -1.0;
  i++;
  (*vertexCoords)[numDims*i+0] = -phi;
  (*vertexCoords)[numDims*i+1] =  0.0;
  (*vertexCoords)[numDims*i+2] =  1.0;
  i++;
  (*vertexCoords)[numDims*i+0] = -phi;
  (*vertexCoords)[numDims*i+1] =  0.0;
  (*vertexCoords)[numDims*i+2] = -1.0;

  /* Apply dither:
  */
  for(i_ul=0;i_ul<*numPoints;i_ul++){
    for(j=0;j<numDims;j++)
      (*vertexCoords)[numDims*i_ul+j] += ditherFrac*(2.0*gsl_rng_uniform(randGen) - 1.0);
  }

  *cells = malloc(sizeof(**cells)*(*numCells));

  i = 0;
  (*cells)[i].id = (unsigned long)i;
  (*cells)[i].vertx[0] = 0;
  (*cells)[i].vertx[1] = 1;
  (*cells)[i].vertx[2] = 3;
  (*cells)[i].vertx[3] = 11;
  (*cells)[i].neigh[0] = NULL; /* An entry ==NULL flags an external face. */
  (*cells)[i].neigh[1] = &(*cells)[7];
  (*cells)[i].neigh[2] = &(*cells)[4];
  (*cells)[i].neigh[3] = &(*cells)[1];

  i = 1;
  (*cells)[i].id = (unsigned long)i;
  (*cells)[i].vertx[0] = 0;
  (*cells)[i].vertx[1] = 1;
  (*cells)[i].vertx[2] = 9;
  (*cells)[i].vertx[3] = 3;
  (*cells)[i].neigh[0] = NULL; /* An entry ==NULL flags an external face. */
  (*cells)[i].neigh[1] = &(*cells)[9];
  (*cells)[i].neigh[2] = &(*cells)[0];
  (*cells)[i].neigh[3] = &(*cells)[2];

  i = 2;
  (*cells)[i].id = (unsigned long)i;
  (*cells)[i].vertx[0] = 0;
  (*cells)[i].vertx[1] = 1;
  (*cells)[i].vertx[2] = 5;
  (*cells)[i].vertx[3] = 9;
  (*cells)[i].neigh[0] = NULL; /* An entry ==NULL flags an external face. */
  (*cells)[i].neigh[1] = &(*cells)[15];
  (*cells)[i].neigh[2] = &(*cells)[1];
  (*cells)[i].neigh[3] = &(*cells)[3];

  i = 3;
  (*cells)[i].id = (unsigned long)i;
  (*cells)[i].vertx[0] = 0;
  (*cells)[i].vertx[1] = 1;
  (*cells)[i].vertx[2] = 7;
  (*cells)[i].vertx[3] = 5;
  (*cells)[i].neigh[0] = NULL; /* An entry ==NULL flags an external face. */
  (*cells)[i].neigh[1] = &(*cells)[13];
  (*cells)[i].neigh[2] = &(*cells)[2];
  (*cells)[i].neigh[3] = &(*cells)[4];

  i = 4;
  (*cells)[i].id = (unsigned long)i;
  (*cells)[i].vertx[0] = 0;
  (*cells)[i].vertx[1] = 1;
  (*cells)[i].vertx[2] = 11;
  (*cells)[i].vertx[3] = 7;
  (*cells)[i].neigh[0] = NULL; /* An entry ==NULL flags an external face. */
  (*cells)[i].neigh[1] = &(*cells)[5];
  (*cells)[i].neigh[2] = &(*cells)[3];
  (*cells)[i].neigh[3] = &(*cells)[0];

  i = 5;
  (*cells)[i].id = (unsigned long)i;
  (*cells)[i].vertx[0] = 0;
  (*cells)[i].vertx[1] = 7;
  (*cells)[i].vertx[2] = 11;
  (*cells)[i].vertx[3] = 12;
  (*cells)[i].neigh[0] = NULL; /* An entry ==NULL flags an external face. */
  (*cells)[i].neigh[1] = &(*cells)[6];
  (*cells)[i].neigh[2] = &(*cells)[12];
  (*cells)[i].neigh[3] = &(*cells)[4];

  i = 6;
  (*cells)[i].id = (unsigned long)i;
  (*cells)[i].vertx[0] = 0;
  (*cells)[i].vertx[1] = 8;
  (*cells)[i].vertx[2] = 12;
  (*cells)[i].vertx[3] = 11;
  (*cells)[i].neigh[0] = NULL; /* An entry ==NULL flags an external face. */
  (*cells)[i].neigh[1] = &(*cells)[5];
  (*cells)[i].neigh[2] = &(*cells)[7];
  (*cells)[i].neigh[3] = &(*cells)[19];

  i = 7;
  (*cells)[i].id = (unsigned long)i;
  (*cells)[i].vertx[0] = 0;
  (*cells)[i].vertx[1] = 3;
  (*cells)[i].vertx[2] = 8;
  (*cells)[i].vertx[3] = 11;
  (*cells)[i].neigh[0] = NULL; /* An entry ==NULL flags an external face. */
  (*cells)[i].neigh[1] = &(*cells)[6];
  (*cells)[i].neigh[2] = &(*cells)[0];
  (*cells)[i].neigh[3] = &(*cells)[8];

  i = 8;
  (*cells)[i].id = (unsigned long)i;
  (*cells)[i].vertx[0] = 0;
  (*cells)[i].vertx[1] = 3;
  (*cells)[i].vertx[2] = 6;
  (*cells)[i].vertx[3] = 8;
  (*cells)[i].neigh[0] = NULL; /* An entry ==NULL flags an external face. */
  (*cells)[i].neigh[1] = &(*cells)[18];
  (*cells)[i].neigh[2] = &(*cells)[7];
  (*cells)[i].neigh[3] = &(*cells)[9];

  i = 9;
  (*cells)[i].id = (unsigned long)i;
  (*cells)[i].vertx[0] = 0;
  (*cells)[i].vertx[1] = 3;
  (*cells)[i].vertx[2] = 9;
  (*cells)[i].vertx[3] = 6;
  (*cells)[i].neigh[0] = NULL; /* An entry ==NULL flags an external face. */
  (*cells)[i].neigh[1] = &(*cells)[16];
  (*cells)[i].neigh[2] = &(*cells)[8];
  (*cells)[i].neigh[3] = &(*cells)[1];

  i = 10;
  (*cells)[i].id = (unsigned long)i;
  (*cells)[i].vertx[0] = 0;
  (*cells)[i].vertx[1] = 2;
  (*cells)[i].vertx[2] = 4;
  (*cells)[i].vertx[3] = 10;
  (*cells)[i].neigh[0] = NULL; /* An entry ==NULL flags an external face. */
  (*cells)[i].neigh[1] = &(*cells)[17];
  (*cells)[i].neigh[2] = &(*cells)[14];
  (*cells)[i].neigh[3] = &(*cells)[11];

  i = 11;
  (*cells)[i].id = (unsigned long)i;
  (*cells)[i].vertx[0] = 0;
  (*cells)[i].vertx[1] = 2;
  (*cells)[i].vertx[2] = 12;
  (*cells)[i].vertx[3] = 4;
  (*cells)[i].neigh[0] = NULL; /* An entry ==NULL flags an external face. */
  (*cells)[i].neigh[1] = &(*cells)[19];
  (*cells)[i].neigh[2] = &(*cells)[10];
  (*cells)[i].neigh[3] = &(*cells)[12];

  i = 12;
  (*cells)[i].id = (unsigned long)i;
  (*cells)[i].vertx[0] = 0;
  (*cells)[i].vertx[1] = 2;
  (*cells)[i].vertx[2] = 7;
  (*cells)[i].vertx[3] = 12;
  (*cells)[i].neigh[0] = NULL; /* An entry ==NULL flags an external face. */
  (*cells)[i].neigh[1] = &(*cells)[5];
  (*cells)[i].neigh[2] = &(*cells)[11];
  (*cells)[i].neigh[3] = &(*cells)[13];

  i = 13;
  (*cells)[i].id = (unsigned long)i;
  (*cells)[i].vertx[0] = 0;
  (*cells)[i].vertx[1] = 2;
  (*cells)[i].vertx[2] = 5;
  (*cells)[i].vertx[3] = 7;
  (*cells)[i].neigh[0] = NULL; /* An entry ==NULL flags an external face. */
  (*cells)[i].neigh[1] = &(*cells)[3];
  (*cells)[i].neigh[2] = &(*cells)[12];
  (*cells)[i].neigh[3] = &(*cells)[14];

  i = 14;
  (*cells)[i].id = (unsigned long)i;
  (*cells)[i].vertx[0] = 0;
  (*cells)[i].vertx[1] = 2;
  (*cells)[i].vertx[2] = 10;
  (*cells)[i].vertx[3] = 5;
  (*cells)[i].neigh[0] = NULL; /* An entry ==NULL flags an external face. */
  (*cells)[i].neigh[1] = &(*cells)[15];
  (*cells)[i].neigh[2] = &(*cells)[13];
  (*cells)[i].neigh[3] = &(*cells)[10];

  i = 15;
  (*cells)[i].id = (unsigned long)i;
  (*cells)[i].vertx[0] = 0;
  (*cells)[i].vertx[1] = 5;
  (*cells)[i].vertx[2] = 10;
  (*cells)[i].vertx[3] = 9;
  (*cells)[i].neigh[0] = NULL; /* An entry ==NULL flags an external face. */
  (*cells)[i].neigh[1] = &(*cells)[16];
  (*cells)[i].neigh[2] = &(*cells)[2];
  (*cells)[i].neigh[3] = &(*cells)[14];

  i = 16;
  (*cells)[i].id = (unsigned long)i;
  (*cells)[i].vertx[0] = 0;
  (*cells)[i].vertx[1] = 6;
  (*cells)[i].vertx[2] = 9;
  (*cells)[i].vertx[3] = 10;
  (*cells)[i].neigh[0] = NULL; /* An entry ==NULL flags an external face. */
  (*cells)[i].neigh[1] = &(*cells)[15];
  (*cells)[i].neigh[2] = &(*cells)[17];
  (*cells)[i].neigh[3] = &(*cells)[9];

  i = 17;
  (*cells)[i].id = (unsigned long)i;
  (*cells)[i].vertx[0] = 0;
  (*cells)[i].vertx[1] = 4;
  (*cells)[i].vertx[2] = 6;
  (*cells)[i].vertx[3] = 10;
  (*cells)[i].neigh[0] = NULL; /* An entry ==NULL flags an external face. */
  (*cells)[i].neigh[1] = &(*cells)[16];
  (*cells)[i].neigh[2] = &(*cells)[10];
  (*cells)[i].neigh[3] = &(*cells)[18];

  i = 18;
  (*cells)[i].id = (unsigned long)i;
  (*cells)[i].vertx[0] = 0;
  (*cells)[i].vertx[1] = 4;
  (*cells)[i].vertx[2] = 8;
  (*cells)[i].vertx[3] = 6;
  (*cells)[i].neigh[0] = NULL; /* An entry ==NULL flags an external face. */
  (*cells)[i].neigh[1] = &(*cells)[8];
  (*cells)[i].neigh[2] = &(*cells)[17];
  (*cells)[i].neigh[3] = &(*cells)[19];

  i = 19;
  (*cells)[i].id = (unsigned long)i;
  (*cells)[i].vertx[0] = 0;
  (*cells)[i].vertx[1] = 4;
  (*cells)[i].vertx[2] = 12;
  (*cells)[i].vertx[3] = 8;
  (*cells)[i].neigh[0] = NULL; /* An entry ==NULL flags an external face. */
  (*cells)[i].neigh[1] = &(*cells)[6];
  (*cells)[i].neigh[2] = &(*cells)[18];
  (*cells)[i].neigh[3] = &(*cells)[11];

  for(i=0;i<(int)(*numCells);i++){
    for(j=0;j<numDims;j++)
      (*cells)[i].centre[j] = 0.0;

    for(k=0;k<numDims+1;k++){
      vi = (*cells)[i].vertx[k];
      for(j=0;j<numDims;j++)
        (*cells)[i].centre[j] += (*vertexCoords)[numDims*vi+j];
    }

    for(j=0;j<numDims;j++)
      (*cells)[i].centre[j] *= oneOnNumVertices;
  }

  *numSides = 30;
  *sides = malloc(sizeof(**sides)*2*(*numSides));

  i = 0;
  (*sides)[2*i  ] = 1;
  (*sides)[2*i+1] = 11;

  i = 1;
  (*sides)[2*i  ] = 11;
  (*sides)[2*i+1] = 3;

  i = 2;
  (*sides)[2*i  ] = 3;
  (*sides)[2*i+1] = 1;

  i = 3;
  (*sides)[2*i  ] = 1;
  (*sides)[2*i+1] = 9;

  i = 4;
  (*sides)[2*i  ] = 1;
  (*sides)[2*i+1] = 5;

  i = 5;
  (*sides)[2*i  ] = 1;
  (*sides)[2*i+1] = 7;

  i = 6;
  (*sides)[2*i  ] = 11;
  (*sides)[2*i+1] = 7;

  i = 7;
  (*sides)[2*i  ] = 11;
  (*sides)[2*i+1] = 12;

  i = 8;
  (*sides)[2*i  ] = 11;
  (*sides)[2*i+1] = 8;

  i = 9;
  (*sides)[2*i  ] = 3;
  (*sides)[2*i+1] = 8;

  i = 10;
  (*sides)[2*i  ] = 3;
  (*sides)[2*i+1] = 6;

  i = 11;
  (*sides)[2*i  ] = 3;
  (*sides)[2*i+1] = 9;

  i = 12;
  (*sides)[2*i  ] = 5;
  (*sides)[2*i+1] = 7;

  i = 13;
  (*sides)[2*i  ] = 7;
  (*sides)[2*i+1] = 12;

  i = 14;
  (*sides)[2*i  ] = 12;
  (*sides)[2*i+1] = 8;

  i = 15;
  (*sides)[2*i  ] = 8;
  (*sides)[2*i+1] = 6;

  i = 16;
  (*sides)[2*i  ] = 6;
  (*sides)[2*i+1] = 9;

  i = 17;
  (*sides)[2*i  ] = 9;
  (*sides)[2*i+1] = 5;

  i = 18;
  (*sides)[2*i  ] = 10;
  (*sides)[2*i+1] = 4;

  i = 19;
  (*sides)[2*i  ] = 4;
  (*sides)[2*i+1] = 2;

  i = 20;
  (*sides)[2*i  ] = 2;
  (*sides)[2*i+1] = 10;

  i = 21;
  (*sides)[2*i  ] = 10;
  (*sides)[2*i+1] = 5;

  i = 22;
  (*sides)[2*i  ] = 10;
  (*sides)[2*i+1] = 9;

  i = 23;
  (*sides)[2*i  ] = 10;
  (*sides)[2*i+1] = 6;

  i = 24;
  (*sides)[2*i  ] = 4;
  (*sides)[2*i+1] = 6;

  i = 25;
  (*sides)[2*i  ] = 4;
  (*sides)[2*i+1] = 8;

  i = 26;
  (*sides)[2*i  ] = 4;
  (*sides)[2*i+1] = 12;

  i = 27;
  (*sides)[2*i  ] = 2;
  (*sides)[2*i+1] = 12;

  i = 28;
  (*sides)[2*i  ] = 2;
  (*sides)[2*i+1] = 7;

  i = 29;
  (*sides)[2*i  ] = 2;
  (*sides)[2*i+1] = 5;
}






