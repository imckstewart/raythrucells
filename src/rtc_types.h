#ifndef RTC_TYPES_H
#define RTC_TYPES_H

#include <stdlib.h>
#include <stdio.h>

#include "dims.h" /* Should define the macro N_DIMS which specifies the number of spatial dimensions. Note that this is a _maximum_ value, the actual value 'numDims' which is passed via the function interfaces may be <= N_DIMS. */

#ifndef TRUE
#define TRUE                   (_Bool)1
#endif

#ifndef FALSE
#define FALSE                  (_Bool)0
#endif

/* In the following comments, N is short for N_DIMS, the number of spatial dimensions. */

struct gridPoint {
  unsigned long id;
  int numNeigh;
  struct gridPoint **neigh;
};

/* The following is just a convenience struct to contain a bunch of stuff which are passed very often together to functions in the 2nd-order interpolation routines. It should be initialized in initializeBaryBuf().
*/
typedef struct  {
  unsigned long id;
  unsigned long vertices[2]; /* these are interpreted as in index to any list of vertex quantities supplied. */
} edgeType;

/* NOTE that it is assumed that vertx[i] is opposite the face that abuts with neigh[i] for all i.
*/ 
struct simplex {
  unsigned long id;
  unsigned long vertx[N_DIMS+1]; /* this is interpreted as in index to any list of vertex quantities supplied. */
  unsigned long edges[N_DIMS*(N_DIMS+1)/2]; /* this is interpreted as in index to any list of edgeType structs supplied. */
  double centre[N_DIMS];
  struct simplex *neigh[N_DIMS+1]; /* An entry ==NULL flags an external face. */
};

/* This struct is meant to record all relevant information about the intersection between a ray (defined by a direction unit vector 'dir' and a starting position 'r') and a face of a simplex.
*/
typedef struct {
  int fi;
  /* The index (in the range {0...N}) of the face (and thus of the opposite vertex, i.e. the one 'missing' from the bary[] list of this face). This index should key to members .vertx and .neigh of struct simplex.
  */
  int orientation;
  /* >0 means the ray exits, <0 means it enters, ==0 means the face is parallel to ray.
  */
  double bary[N_DIMS], dist, collPar;
  /* These are the barycentric coordinates of the point of intersection between the ray and the face 'fi' of the simplex. The quantity 'dist' is defined via r_int = r + dist*dir. 'collPar' is a measure of how close r_int lies to any edge of the face. Note that the ordering of the indices of .bary is the same as the ordering for .vertx and .neigh of the simplex it belongs to, only with the vertex corresponding to the intersected face missing.
  */
} intersectType;

typedef struct {
  int nCellsInChain,nCellsMallocd;
  intersectType entryIntcpt,*exitIntcpts;
  unsigned long *cellIds;
  _Bool *flags;
} cellChainType;

typedef struct {
  double r[N_DIMS][N_DIMS], simplexCentre[N_DIMS];
  /* 'r' is a list of the the N vertices of the face, each of which has N cartesian components. 'simplexCentre' is a convenience pointer which gives the location of the geometric centre of the simplex. */
} faceType;

typedef struct {
//  double axes[N_DIMS-1][N_DIMS], r[N_DIMS][N_DIMS-1], origin[N_DIMS];
  double axes[N_DIMS][N_DIMS], r[N_DIMS][N_DIMS-1], origin[N_DIMS];
  /* 'r' expresses the location of the N vertices of a simplicial polytope face in N-space, in terms of components along the N-1 orthogonal axes in the sub-plane of the face. Thus you should malloc r as r[N][N-1]. */
} facePlusBasisType;

typedef struct {
  int cellAlongRayI;
  _Bool pixelIsInCells;
} rasterType;

typedef struct {
  int numPixels;
  double origin, delta;
} axisType;

typedef struct {
  int numDims,numVertices,numEdges,numElements,(*edgeVertexIndices)[2];
  double *vertexValues,*edgeValues,*vertxShapeValues,*edgeShapeValues;
} baryBuffType;



struct gridPoint	init_gridPoint(void);
void			free_gridPoints(struct gridPoint *gps, const unsigned long numPoints);
edgeType		init_edge(void);
struct simplex		init_simplex(void);
cellChainType		init_cellChain(const int nCellsMallocd);
cellChainType*		realloc_cellChain(cellChainType *chain, const int nCellsMallocd);
void			free_cellChain(cellChainType *chain);
intersectType		init_intersect(void);
faceType		init_face(void);
facePlusBasisType	init_facePlusBasis(void);
void			init_baryBuff(const int numDims, const int numElements, baryBuffType *baryBuff);
void			free_baryBuff(baryBuffType *baryBuff);


#endif /* RTC_TYPES_H */

