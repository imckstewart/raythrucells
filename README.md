Project raythrucells
====================

Suppose you have a connected, convex set of simplicial cells in a space with N dimensions. What exactly does this mean? Well, a simplex (to borrow from the Wikipedia article) is a generalization of the notion of a triangle or tetrahedron to arbitrary dimensions. A simplex in N dimensions will have N+1 vertices and N+1 faces, each face of which will be a simplex in N-1 dimensions. 'Connected' means that you can travel between any two simplices in the set without ever leaving a simplex. 'Convex' means that the convex hull (google it) of the set of points of all the simplices contains no space that is not part of a simplex. In simple language: there are no holes or bends in the set of simplices.

So we have this set, and we want to pass a ray through it: i.e we want to know which cells the ray intersects, plus their order; and we want as a bonus to know the spatial locations of the points at which the ray enters or exits each cell. This might be to solve for some property along the ray (radiative transfer equations perhaps). The code in the present project allows you to do that. It is not such a trivial business. Any kind of computational geometry runs up against the same problems: that certain conditions may be absolutely determinable gemoetrically (e.g. whether a point lies within a polytope) but, due to rounding errors, not determinable with complete accuracy in computing. The present algorithm chases up all the borderline cases until it reaches the best resolution of the problem.

Any C program which wishes to use this code should look broadly like the following:

```
#include "raythrucells.h"

int main() {
  int status=0;
  double x[N_DIMS],dir[N_DIMS],*vertexCoords=NULL;
  struct simplex *cells=NULL;
  unsigned long numCells;
  const double epsilon=1.0e-6;
  intersectType entryIntcpt,*cellExitIntcpts=NULL;
  unsigned long *chainOfCellIds=NULL;
  int lenChainPtrs;

  < Construct the mesh of cells, which will involve mallocing 'cells' and 'vertexCoords' and setting their elements appropriately. >

  status = followRayThroughCells(N_DIMS, x, dir, cells, numCells, epsilon\
    , NULL, vertexCoords, &entryIntcpt, &chainOfCellIds, &cellExitIntcpts\
    , &lenChainPtrs);

  if(status!=0){
    printf("followRayThroughCells status return %d\n", status);
exit(1);
  } // Note that a status value of RTC_ERR_NO_ENTRIES indicates that the ray did not intersect with any of the cells.

  < Make use of the information. >

  free(chainOfCellIds);
  free(cellExitIntcpts);
  free(vertexCoords);
  free(cells);

  return 0;
}
```

