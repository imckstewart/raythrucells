#include "raythrucells.h"

#ifdef TEST
#include "../test/test.h"
#endif

/*....................................................................*/
faceType
_extractFace(const int numDims, double *vertexCoords, struct simplex *cells\
  , const unsigned long dci, const int fi){
  /*
Given a simplex cells[dci] and the face index (in the range {0...numDims}) fi, this returns the desired information about that face. Note that the ordering of the elements of face.r[] is the same as the ordering of the vertices of the simplex, cells[dci].vertx[]; just the vertex fi is omitted.

Note that the element 'centre' of the faceType struct is mean to contain the spatial coordinates of the centre of the simplex, not of the face. This is designed to facilitate orientation of the face and thus to help determine whether rays which cross it are entering or exiting the simplex.
 */

  const int numFaces=numDims+1;
  int vi,vvi,di;
  faceType face;
  unsigned long gi;

  vvi = 0;
  for(vi=0;vi<numFaces;vi++){
    if(vi!=fi){
      gi = cells[dci].vertx[vi];
      for(di=0;di<numDims;di++){
        face.r[vvi][di] = vertexCoords[numDims*gi+di];
      }
      vvi++;
    }
  }

  for(di=0;di<numDims;di++)
    face.simplexCentre[di] = cells[dci].centre[di];

return face;
}

/*....................................................................*/
int
_getNewEntryFaceI(const int numDims, const unsigned long dci, const struct simplex *newCell){
  /* Finds the index of the old cell in the face list of the new cell. */

  const int numFaces=numDims+1;
  _Bool matchFound = FALSE;
  int ffi = 0, newEntryFaceI;

  while(ffi<numFaces && !matchFound){
    if(newCell->neigh[ffi]!=NULL && newCell->neigh[ffi]->id==dci){
      matchFound = TRUE;
      newEntryFaceI = ffi;
    }
    ffi++;
  }

  /* Sanity check: */
  if(!matchFound)
    rtcError(RTC_ERR_OLD_NOT_FOUND, "Cannot find old cell ID in new cell data.");

return newEntryFaceI;
}

/*....................................................................*/
facePlusBasisType
_calcFaceInNMinus1(const int numDims, faceType *face){
  /*
Each of the faces of a simplex in N spatial dimensions is itself a simplex in N-1 dimensions. Each face can therefore be represented via the coordinates of its vertices expressed in an (N-1)-dimensional frame oriented so as to be parallel to the face. The function of the present routine is to perform this decomposition and to return both the (N-1)-dimensional frame (specified via the coordinates in N dimensions of its N-1 basis vectors) and the coordinates in it of each of the N face vertices.
  */

  facePlusBasisType facePlusBasis;
  double edgesMinusOrigin[numDims-1][numDims],dotValue;
  int di,vi,ddi;

  /* The first part of the routine is finding the components in N-space of the N-1 orthonormal axes parallel to the face.
  */
  for(di=0;di<numDims;di++)
    facePlusBasis.origin[di] = face->r[0][di];

  for(vi=0;vi<numDims-1;vi++){
    for(di=0;di<numDims;di++)
      edgesMinusOrigin[vi][di] = face->r[vi+1][di] - facePlusBasis.origin[di];
  }

  gramSchmidt(numDims, numDims-1, edgesMinusOrigin, facePlusBasis.axes);

  /* Now we calculate the coords of each vertex in the N-1 system:
  */
  vi = 0;
  for(ddi=0;ddi<numDims-1;ddi++)
    facePlusBasis.r[vi][ddi] = 0.0;

  for(vi=1;vi<numDims;vi++){
    for(ddi=0;ddi<numDims-1;ddi++){
      dotValue = calcDotProduct(numDims, edgesMinusOrigin[vi-1], facePlusBasis.axes[ddi]);
      facePlusBasis.r[vi][ddi] = dotValue;
    }
  }

  return facePlusBasis;
}

/*....................................................................*/
intersectType
_intersectLineWithFace(const int numDims, double *x, double *dir, faceType *face\
  , const double epsilon){
  /*
This function calculates the intersection between a line and the face of a simplex oriented in that space. Obviously the number of dimensions of the space must be >=2. The intersection point may be expressed as

	px_ = x_ + a*dir_

where px_, x_ and dir_ are vectors and 'a' is a scalar. The routine returns the following information:
	- The value of 'a' (returned in field 'dist' of the returned intersectType).
	- The so-called barycentric coordinates (BC) of px_ in the face.

The scalar 'a' is found as follows. We need the additional point y_ which can be any of the face vertices. We state the vector identity

	a*dir_ = (y_ - x_) + (px_ - y_).

Clearly (px_ - y_) is parallel to the face. If we take the scalar product of both sides with the vector n_ which is normal to the face we arrive at

	a*dir_.n_ = (y_ - x_).n_ + (px_ - y_).n_.

	          = (y_ - x_).n_ + 0.
Thus
	     (y_ - x_).n_
	a = --------------.
	       dir_.n_

Notes:
	* This routine works best when the sides of the face are not too disparate in size.

	* There is, of course, no guarantee that the line actually intersects the face, even if the line and the face are non-parallel. There are also borderline cases, i.e. where the line passes close to a vertex, in which an exact calculation would show that the intersection occurs (or doesn't occur), but the imprecise computed value claims that it doesn't (or does). Intersection may be judged via the values of the barycentric coordinates (BC): if the BC all fall in the interval [0,1], the line intersects the face; if one of the BC is negative, it doesn't.

        * The routine does not fill in the value of intcpt.fi because it does not have access to vertex indices.
  */
  const double oneOnEpsilon=1.0/epsilon;
  double vs[numDims-1][numDims],norm[numDims],normDotDx,numerator,pxInFace[numDims-1];
  int i,j,k,di,ci,ri,vi,ciOfMax,ciOfMin;
  double testSumForCW=0.0,maxSingularValue,singularValue;
  facePlusBasisType facePlusBasis;
  char errStr[RTC_MSG_STR_LEN];
  intersectType intcpt;

  for(vi=0;vi<numDims-1;vi++){
    for(di=0;di<numDims;di++)
      vs[vi][di] = (*face).r[vi+1][di]-(*face).r[0][di];
  }

  /* First calculate a normal vector to the face (note that it doesn't need to be of length==1).
  */
  if(numDims==2){
    norm[0] = -vs[0][1];
    norm[1] =  vs[0][0];

  }else if(numDims==3){
    /* Calculate norm via cross product. */
    for(di=0;di<numDims;di++){
      j = (di+1)%numDims;
      k = (di+2)%numDims;
      norm[di] = vs[0][j]*vs[1][k] - vs[0][k]*vs[1][j];
    }

  }else{ /* Assume numDims>3 */
    /* Calculate norm via SVD. */
    int status=0;
    gsl_matrix *matrix = gsl_matrix_alloc(numDims-1, numDims);
    gsl_matrix *svv    = gsl_matrix_alloc(numDims,   numDims);
    gsl_vector *svs    = gsl_vector_alloc(numDims);
    gsl_vector *work   = gsl_vector_alloc(numDims);

    for(ci=0;ci<numDims;ci++){
      for(ri=0;ri<numDims-1;ri++)
        gsl_matrix_set(matrix, ri, ci, vs[ri][ci]);
    }

    status = gsl_linalg_SV_decomp(matrix, svv, svs, work);
    if(status){
      snprintf(errStr, RTC_MSG_STR_LEN, "SVD decomposition failed (GSL error %d).", status);
      rtcError(RTC_ERR_SVD_FAIL, errStr);
    }

    /*
Since we have N-1 equations in N unknowns, we would expect at least one of the N elements of svs to be zero (within rounding error). The column of svv which corresponds to this value should then be the normal vector which we require. We'll just check however that not more than 1 value of svs is zero.

The GSL doco says that SV_decomp returns sorted svs values, but I prefer not to rely on that.
    */

    ci = 0;
    ciOfMax = ci;
    maxSingularValue = gsl_vector_get(svs,ci);
    for(ci=1;ci<numDims;ci++){
      singularValue = gsl_vector_get(svs,ci);
      if(singularValue>maxSingularValue){
        ciOfMax = ci;
        maxSingularValue = singularValue;
      }
    }

    ciOfMin = -1; /* Impossible default. */
    for(ci=0;ci<numDims;ci++){
      if(ci==ciOfMax) continue;

      singularValue = gsl_vector_get(svs,ci);
      if(singularValue*oneOnEpsilon<maxSingularValue){
        if(ciOfMin>=0){
          /* This is an error because it indicates that >1 singular values are 'small'. */
          snprintf(errStr, RTC_MSG_STR_LEN, "Simplex face does not span an N-1 subspace.");
          rtcError(RTC_ERR_NON_SPAN, errStr);
        }

        ciOfMin = ci;
      }
    }

    for(di=0;di<numDims;di++)
      norm[di] = gsl_matrix_get(svv,di,ciOfMin);

    gsl_vector_free(work);
    gsl_vector_free(svs);
    gsl_matrix_free(svv);
    gsl_matrix_free(matrix);
  }

  /* Since we don't know a priori whether the vertices of the face are listed CW or ACW seen from inside the simplex, we will work this out by dotting the normal vector with a vector from the centre of the simplex to vertex 0 of the face. (A simplex is always convex so this ought to work.)
  */
  testSumForCW = 0.0;
  for(di=0;di<numDims;di++)
    testSumForCW += norm[di]*((*face).r[0][di] - (*face).simplexCentre[di]);

  if(testSumForCW<0.0){
    for(di=0;di<numDims;di++)
      norm[di] *= -1.0;
  }

  /* Calculate the scalar (or dot) product between norm and dir.
  */
  normDotDx = 0.0;
  for(di=0;di<numDims;di++)
    normDotDx += norm[di]*dir[di];

  if(normDotDx>0.0){ /* it is an exit face. */
    intcpt.orientation = 1;
  }else if(normDotDx<0.0){ /* it is an entry face. */
    intcpt.orientation = -1;
  }else{ /* normDotDx==0.0, i.e. line and face are parallel. */
    intcpt.orientation = 0;
    for(di=0;di<numDims;di++)
      intcpt.bary[di] = 0.0;
    intcpt.dist    = 0.0;
    intcpt.collPar = 0.0;

    return intcpt;
  }

  /* If we've got to here, we can be sure the line and the face are not parallel, and that we therefore expect meaningful results for the calculation of 'a'.
  */
  numerator = 0.0;
  for(di=0;di<numDims;di++)
    numerator += norm[di]*((*face).r[0][di] - x[di]); /* n_.(y_ - x_) */

  intcpt.dist = numerator/normDotDx;

  /* In order to calculate the barycentric coordinates, we need to set up a N-1 coordinate basis in the plane of the face.
  */
  facePlusBasis = _calcFaceInNMinus1(numDims, face);

  /* Now we want to express the intersection point in these coordinates:
  */
  for(i=0;i<numDims-1;i++){
    pxInFace[i] = 0.0;
    for(di=0;di<numDims;di++)
      pxInFace[i] += (x[di] + intcpt.dist*dir[di] - facePlusBasis.origin[di])*facePlusBasis.axes[i][di];
  }

  /*
The barycentric coordinates x_ = {L_1,L_2,...,L_{N-1}} are given by

	T x_ = b_

where T is an (N-1)*(N-1) matrix with entries

	T_{i,j} = facePlusBasis.r[j+1][i] - facePlusBasis.r[0][i]

and
	b_i = pxInFace[i] - facePlusBasis.r[0][i].

The final BC L_0 is given by

	         _N-1
	         \
	L_0 = 1 - >    L_i.
	         /_i=1
  */

  calcBaryCoords(numDims-1, facePlusBasis.r, pxInFace, intcpt.bary);

  /* Finally, calculate the 'collision parameter':
  */
  di = 0;
  if(intcpt.bary[di] < 0.5)
    intcpt.collPar = intcpt.bary[di];
  else
    intcpt.collPar = 1.0 - intcpt.bary[di];

  for(di=1;di<numDims;di++){
    if(intcpt.bary[di] < 0.5){
      if(intcpt.bary[di] < intcpt.collPar)
        intcpt.collPar = intcpt.bary[di];
    }else{ /* intcpt.bary[di]>=0.5 */
      if(1.0 - intcpt.bary[di] < intcpt.collPar)
        intcpt.collPar = 1.0 - intcpt.bary[di];
    }
  }

  return intcpt;
}

/*....................................................................*/
_Bool
_followGoodChain(const int numDims, double *x, double *dir\
  , struct simplex *cells, unsigned long dci, int cellEntryFaceI\
  , const double epsilon, faceType **facePtrs[numDims+1]\
  , double *vertexCoords, _Bool **cellVisited, cellChainType *cellChain\
  , intersectType marginalExitIntcpts[numDims+1], int *numMarginalExits){
  /*
We're following a ray through a set of cells. The end aim is to produce several lists (stored in the struct 'cellChain') pertaining to the N cells the ray passes through. We want a list of the N cell IDs and a list of the N+1 intercepts between the ray and each of the cell faces it encounters.

The present function follows the ray through a succession of cells only so long as the decision about which next cell is entered remains unambiguous - or, until the far edge of the set is encountered. (There's at present no provision for non-convex sets.)

The function returns TRUE if the chain terminated at the edge of the set of cells, following only good interceptions.

Input arguments:
================
	- numDims: the number of spatial dimensions.

	- x: the nominal origin of the ray.

	- dir: a unit vector giving the direction of propagation of the ray.

	- cells: a list of all the cells in the connected set.

	- dci: at entry to the function, this is the ID of the first cell the function will process. The function changes its value, but since it is a variable with only local scope, the changed value is not returned.

	- cellEntryFaceI: at entry to the function, this gives the index (in the range [0,numDims]) of the entry face in the list of faces of the entry cell. The function changes its value, but since it is a variable with only local scope, the changed value is not returned.

	- epsilon: a small value.

	- facePtrs: contains information about the spatial locations of all vertices in the set of cells. Either this or vertexCoords can be NULL. See explanation in followRayThroughCells().

	- vertexCoords: contains information about the spatial locations of all vertices in the set of cells. Either this or facePtrs can be NULL. See explanation in followRayThroughCells().

Input/output arguments:
=======================
	- cellVisited: self-explanatory.

	- cellChain: this contains information about the cells which are traversed by the ray.

Output arguments:
=================
	- marginalExitIntcpts: this is a list of intercepts between the ray and all of the faces of the last cell processed, with the proviso that the ray must be exiting the cell, and that the ID of the new cell is ambiguous (as can happen in numerical geometry).

	- numMarginalExits: size of marginalExitIntcpts.
  */

  const int numFaces=numDims+1;
  int numGoodExits,fi,goodExitFis[numFaces],marginalExitFis[numFaces],exitFi;
  faceType face;
  intersectType cellIntcpt[numFaces];

  while(TRUE){ /* We'll break out of it (via function returns) either if the chain terminates at the edge of the set of cells or as soon as we have no good cell exits and no single marginal exit. */
    (*cellVisited)[dci] = TRUE;

    /* If there is not enough room in chainOfCellIds and cellExitIntcpts, realloc them to new value of lenChainPtrs. */
    if(cellChain->nCellsInChain >= cellChain->nCellsMallocd)
      cellChain = realloc_cellChain(cellChain, cellChain->nCellsMallocd + RTC_BUFFER_SIZE);

    /* Store the current cell ID (we leave storing the exit face for later, when we know what it is). */
    cellChain->cellIds[cellChain->nCellsInChain] = dci;

    /* Calculate numbers of good and marginal exits.
    */
    numGoodExits = 0;
    *numMarginalExits = 0;

    for(fi=0;fi<numFaces;fi++){
      if(fi!=cellEntryFaceI && (cells[dci].neigh[fi]==NULL || !(*cellVisited)[cells[dci].neigh[fi]->id])){
        /* Store points for this face: */
        if(facePtrs==NULL){
          face = _extractFace(numDims, vertexCoords, cells, dci, fi);
        }else{
          face = (*facePtrs)[dci][fi];
        }

        /* Now calculate the intercept: */
        cellIntcpt[fi] = _intersectLineWithFace(numDims, x, dir, &face, epsilon);
        cellIntcpt[fi].fi = fi; /* Ultimately we need this so we can relate the bary coords for the face back to the Delaunay cell. */

        if(cellIntcpt[fi].orientation>0){ /* it is an exit face. */
          if(cellIntcpt[fi].collPar-epsilon>0.0){
            goodExitFis[numGoodExits] = fi;
            numGoodExits++;
          }else if (cellIntcpt[fi].collPar+epsilon>0.0){
            marginalExitFis[*numMarginalExits] = fi;
            marginalExitIntcpts[*numMarginalExits] = cellIntcpt[fi];
            (*numMarginalExits)++;
          }
        }
      }
    }

    if(numGoodExits>1)
      rtcError(RTC_ERR_BUG, "Some sort of bug: more than 1 firm candidate found for ray exit from cell.");

    if(numGoodExits<1 && (*numMarginalExits)!=1)
return FALSE;

    /* If we have reached here then we only have a single exit from the present cell (rated either 'good' or 'marginal').
    */
    if(numGoodExits==1)
      exitFi = goodExitFis[0];
    else /* (*numMarginalExits)==1 */
      exitFi = marginalExitFis[0];

    /* Store the exit face details: */
    cellChain->exitIntcpts[cellChain->nCellsInChain] = cellIntcpt[exitFi];

    cellChain->nCellsInChain++;

    if(cells[dci].neigh[exitFi]==NULL){ /* Signals that we have reached the edge of the model. */
      /* Realloc the ptrs to their final sizes:
      */
      cellChain = realloc_cellChain(cellChain, cellChain->nCellsInChain);
return TRUE;
    }

    cellEntryFaceI = _getNewEntryFaceI(numDims, dci, cells[dci].neigh[exitFi]);
    dci = cells[dci].neigh[exitFi]->id;
  };
}

/*....................................................................*/
int
_buildRayCellChain(const int numDims, double *x, double *dir\
  , struct simplex *cells, unsigned long dci\
  , int entryFaceI, int levelI, const double epsilon\
  , faceType **facePtrs[numDims+1], double *vertexCoords\
  , _Bool **cellVisited, cellChainType *cellChain){
  /*
This function is designed to follow a ray (defined by a starting locus 'x' and a direction vector 'dir') through a convex connected set of cells (assumed simplicial). The function returns an integer status value directly, and two lists (plus their common length) via the argument interface: chainOfCellIds and cellExitIntcpts. Taken together, these lists define a chain of cells traversed by the ray.

The task of determining which cells are traversed by the ray is simple in principle, but complications arise in computational practice due to the finite precision of floating-point calculations. Where the ray crosses a cell face near to one of its edges, numerical calculation of the 'impact parameter' may return an answer which is erroneous either way: i.e., a ray which just misses a face may be reported as hitting it, and vice versa. To deal with this, a distinction is made between impacts which are (i) 'good', that is far from any face edge; (ii) 'bad', i.e. clearly missing the face; and (iii) 'marginal', i.e. closer to the edge than some preset cutoff which is represented in the argument list by the number 'epsilon'. Note that a tally is kept of those cells which have already been tested for the current ray, and any face which abuts a neighbouring cell which has been visited already will be flagged as 'bad'.

The function therefore looks at all the exit faces of the cell and sorts them into these three categories. How it proceeds then depends on the relative numbers of each, as described below.

	- If there is more than 1 good exit face, an exception is generated. This is a sign that epsilon has been chosen with too small a value.

	- If there is just 1 good exit face, the marginal ones are ignored. The current cell data are appended to the chain and the cell abutting the exit face becomes the new working cell.

	- If there are no good exit faces, we look at the marginal ones. If there is only 1 of these, we loop as before. If there are more than 1, the function is called recursively for each cell on the far side of an exit face.

Thus there are two alternate modes of operation within the function: a straightforward loop along the cells in a single chain, which will continue so long as there is only a single exit face of type either 'good' or 'marginal'; and a recursive launch of the function at a fork in the chain into each of its several possible branches.

The function terminates under the following conditions:
	- It detects that the edge of the model is reached (returns success).
	- There are no exit faces (returns failure).
	- There are no good exit faces and either
		* all of the recursive calls to marginal faces have been unsuccessful (returns failure), or
		* one of these has been successful (returns success).

At a successful termination, therefore, details of all the cells to the edge of the model are correctly stored in chainOfCellIds and cellExitIntcpts, and the number of these cells is returned in lenChainPtrs.

The arguments are the same as those of _followGoodChain().

***** Note that it is assumed here that a mis-indentification of the actual cell traversed by a ray in marginal cases will not ultimately matter to the external calling routine, provided that the ray-face intersection points are sufficiently close to the boundary between the cells in question. This is only reasonable if whatever function or property is being sampled by the ray does not vary in a stepwise manner at any cell boundary. *****
  */

  const int numFaces=numDims+1;
  _Bool chainEndedOk;
  int numMarginalExits,exitFi,i,status,newEntryFaceI;
  intersectType marginalExitIntcpts[numFaces];
  unsigned long lastGoodCellId;

  chainEndedOk = _followGoodChain(numDims, x, dir, cells, dci\
    , entryFaceI, epsilon, facePtrs, vertexCoords, cellVisited, cellChain\
    , marginalExitIntcpts, &numMarginalExits);

  if (chainEndedOk)
return 0;

  /* If we got to here, we have run out of good (or at least single) exit-face options, let's try the marginal ones. */

  if(numMarginalExits<1)
return RTC_ERR_BAD_CHAIN; /* Unsuccessful end of this chain. */

  /* If we have got to this point, we must have numMarginalExits>1; thus we have a fork in the chain, and must explore each branch. We recurse here because a recursive scheme is the best way to do that.
  */
  lastGoodCellId = cellChain->cellIds[cellChain->nCellsInChain];
  for(i=0;i<numMarginalExits;i++){
    cellChain->exitIntcpts[cellChain->nCellsInChain] = marginalExitIntcpts[i];
    exitFi = marginalExitIntcpts[i].fi;

    if(cells[lastGoodCellId].neigh[exitFi]==NULL){ /* Signals that we have reached the edge of the model. */
      /* Realloc the ptrs to their final sizes: */
      cellChain->nCellsInChain++;
      cellChain = realloc_cellChain(cellChain, cellChain->nCellsInChain);
      status = 0;

    }else{
      newEntryFaceI = _getNewEntryFaceI(numDims, lastGoodCellId, cells[lastGoodCellId].neigh[exitFi]);

      /* Now we dive into the branch: */
      status = _buildRayCellChain(numDims, x, dir, cells\
        , cells[lastGoodCellId].neigh[exitFi]->id, newEntryFaceI, levelI+1\
        , epsilon, facePtrs, vertexCoords, cellVisited, cellChain);
    }

    if(status==0) break;
  }

return status;
}

/*....................................................................*/
int
followRayThroughCells(const int numDims, double *x, double *dir\
  , struct simplex *cells, const unsigned long numCells, const double epsilon\
  , faceType **facePtrs[numDims+1], double *vertexCoords\
  , cellChainType *cellChain){
  /*
The present function follows a ray through a connected, convex set of cells (assumed to be simplices) and returns information about the chain of cells it passes through. If the ray is found to pass through 1 or more cells, the function returns 0, indicating success; if not, it returns a non-zero value. The chain description consists of three pieces of information: (i) intercept information for the entry face of the first cell encountered; (ii) the IDs of the cells in the chain; (iii) intercept information for the exit face of the ith cell.

The function arguments:
=======================
    Inputs:
    -------
	numDims: the number of spatial dimensions. Note that this may be <= N_DIM.

	x, dir: these are both assumed to be vectors of length numDims. They define the ray, 'x' being its starting point and 'dir' its direction.

	cells: a list of numCells entries which defines the mesh of connected cells. The coordinates of the vertices of each cell can be obtained in two different ways: (i) by finding the matching entry in the list facePtrs; (ii) from the list vertexCoords, indexed by the element 'vertx' of struct simplex.

	numCells: as it says.

	epsilon: a small value. 1e-6??

	facePtrs: a list, with one entry per cell, which stores information (including the location of the face vertices) about the numDims+1 faces of that cell. The ith face of any cell is assumed to abut onto the cell pointed to by the ith entry of element 'neigh' of the cell struct. It may be set to NULL, but if so, the argument 'vertexCoords' may not be NULL.

	vertexCoords: a vector of size numDims*numPoints, where numPoints is not supplied or used here, but it is assumed that the indices stored in the 'vertx' element of the cells never equals or exceeds it. If 'facePtrs' is supplied (i.e.non-NULL), then 'vertexCoords' is not used and may be left at NULL.

    Outputs:
    --------
	cellChain: contains the following fields of present interest:
		entryIntcpt: this gives the location where the ray (extended backward beyond its given starting point if necessary) first enters the mesh of cells. ***NOTE*** that since it is a scalar, non-convex meshes, for which several points of entry are possible, will not be handled correctly.

		cellIds: this is a list, in order, of the cells traversed by the ray from its first entry to its exit from the mesh of cells.

		exitIntcpts: a list of the locations where the ray exits each cell. Together with 'entryIntcpt', this enables the precise path of the ray through each cell in 'chainOfCellIds' to be known.

		nCellsInChain: the length of the returned lists 'chainOfCellIds' and 'cellExitIntcpts'.


Notes:
======
The pointer *cellChain should be freed (via the function free_cellChain) after the present function is called.

The argument facePtrs may be set to NULL, in which case the function will construct each face from the list of cells etc as it needs it. This saves on memory but takes more time. If the calling routine supplies these values it needs to do something like as follows:

	faceType face,*facePtrs[N_DIMS+1]=malloc(sizeof(*(*facePtrs[N_DIMS+1]))*numCells);
	for(dci=0;dci<numCells;dci++){
	  for(j=0;j<numDims+1;j++){
	    face = __extractFace(numDims, vertexCoords, cells, dci, j);
	    facePtrs[dci][j] = face;
	  }
	}
	status = followRayThroughCells(... &facePtrs, ...);

Note finally that if facePtrs is supplied non-NULL, vertexCoords may be left at NULL. If this is filled, it should be malloc'd as

	vertexCoords = malloc(sizeof(double)*numDims*numPoints);

and filled as

	for(i=0;i<numPoints;i++)
	  for(j=0;j<numDims;j++)
	    vertexCoords[numDims*i+j] = // grid point i, coordinate j

  */

  const int numFaces=numDims+1, maxNumEntryFaces=100;
  int numEntryFaces,fi,entryFis[maxNumEntryFaces],i,status;
  faceType face;
  unsigned long dci,entryDcis[maxNumEntryFaces];
  intersectType intcpt,entryIntcpts[maxNumEntryFaces];
  _Bool *cellVisited=NULL;

  /* Choose a set of starting faces by testing all the 'external' faces of cells which have some. */
  numEntryFaces = 0;
  for(dci=0;dci<numCells;dci++){
    for(fi=0;fi<numFaces;fi++){
      if(cells[dci].neigh[fi]==NULL){ /* means that this face lies on the outside of the model. */
        /* Store points for this face: */
        if(facePtrs==NULL){
          face = _extractFace(numDims, vertexCoords, cells, dci, fi);
        }else{
          face = (*facePtrs)[dci][fi];
        }

        /* Now calculate the intercept: */
        intcpt = _intersectLineWithFace(numDims, x, dir, &face, epsilon);
        intcpt.fi = fi; /* Ultimately we need this so we can relate the bary coords for the face back to the Delaunay cell. */

        if(intcpt.orientation<0){ /* it is an entry face. */
          if(intcpt.collPar+epsilon>0.0){
            if(numEntryFaces>maxNumEntryFaces)
              rtcError(RTC_ERR_TOO_MANY_ENTRY, "Too many entry faces.");

            entryDcis[   numEntryFaces] = dci;
            entryFis[    numEntryFaces] = fi;
            entryIntcpts[numEntryFaces] = intcpt;
            entryIntcpts[numEntryFaces] = intcpt;
            numEntryFaces++;
          }
        }
      }
    }
  }

  if(numEntryFaces<=0){
    *cellChain = init_cellChain(0);

return 0; /* This is ok, it can happen if the ray missed the mesh of cells entirely. */
  }

  *cellChain = init_cellChain(RTC_BUFFER_SIZE);
  cellVisited = malloc(sizeof(*cellVisited)*numCells);
  for(dci=0;dci<numCells;dci++)
    cellVisited[dci] = FALSE;

  for(i=0;i<numEntryFaces;i++){
    status = _buildRayCellChain(numDims, x, dir, cells\
      , entryDcis[i], entryFis[i], 0, epsilon, facePtrs, vertexCoords\
      , &cellVisited, cellChain);

    if(status==0) break;
  }

  if(status==0){ /* means ith entry face returned a good chain. */
    cellChain->entryIntcpt = entryIntcpts[i];
    /* Note that the order of the bary coords, and the value of fi, are with reference to the vertx list of the _entered_ cell. This can't of course be any other way, because this ray enters this first cell from the exterior of the model, where there are no cells. For all the intersectType objects in the list cellExitIntcpts, the bary coords etc are with reference to the exited cell. */

  }else{ /* means none of the possibly >1 entry faces returned a good chain. Return with status value from last bad chain. */
    free_cellChain(cellChain);
    *cellChain = init_cellChain(0);
  }

  free(cellVisited);

return status;
}

