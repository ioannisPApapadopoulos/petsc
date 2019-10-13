static char help[] = "Tests for submesh creation\n\n";

/* TODO
*/

#include <petscdmplex.h>
#include <petscsf.h>
#include <petsc/private/dmpleximpl.h>    /*I      "petscdmplex.h"    I*/
#include <petsc/private/dmlabelimpl.h>   /*I      "petscdmlabel.h"   I*/




static PetscErrorCode _DMPlexMarkSubpointMap_Closure(DM dm, DMLabel filter,
                                                    PetscInt filterValue,
                                                    PetscInt height,
                                                    DMLabel subpointmap,
                                                    PetscBool addOverlap)
{
  PetscErrorCode     ierr;
  IS                 marked;
  const PetscInt    *points;
  PetscInt          *closure = NULL, *pStart, *pEnd;
  PetscInt           p, point, leafpoint, npoints, hStart, hEnd, tdim, d;
  PetscSF            sf;
  PetscInt           nroots, nleaves;
  const PetscInt    *ilocal;
  const PetscSFNode *iremote;
  PetscHMapI         pointmap;
  PetscInt           tempMark, tempMarkGhost, defaultValue;

  PetscFunctionBegin;

  ierr = DMLabelGetDefaultValue(subpointmap, &defaultValue); CHKERRQ(ierr);
  ierr = DMPlexGetHeightStratum(dm, height, &hStart, &hEnd); CHKERRQ(ierr);

  /* Mark non-overlapping points */
  ierr = DMGetDimension(dm, &tdim); CHKERRQ(ierr);
  ierr = DMLabelGetStratumIS(filter, filterValue, &marked); CHKERRQ(ierr);
  tempMark = tdim + 1;
  tempMarkGhost = tempMark + 1;
  ierr = DMGetPointSF(dm, &sf); CHKERRQ(ierr);
  ierr = PetscSFGetGraph(sf, &nroots, &nleaves, &ilocal, &iremote); CHKERRQ(ierr);

  ierr = PetscHMapICreate(&pointmap); CHKERRQ(ierr);
  for (p = 0; p < nleaves; ++p) {
    PetscHMapISet(pointmap, ilocal[p], p);
  }
  if (marked) {
    ierr = ISGetLocalSize(marked, &npoints); CHKERRQ(ierr);
    ierr = ISGetIndices(marked, &points); CHKERRQ(ierr);
    PetscInt closureSize, ci;
    /* Mark/Overwrite closure of core points */
    for (p = 0; p < npoints; ++p) {
      point = points[p];
      if (point < hStart || point >= hEnd) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Filter label marks a point at the incorrect height");
      ierr = PetscHMapIGet(pointmap, point, &leafpoint); CHKERRQ(ierr);
      if (leafpoint < 0) {
        ierr = DMPlexGetTransitiveClosure(dm, point, PETSC_TRUE, &closureSize, &closure); CHKERRQ(ierr);
        for (ci = 0; ci < closureSize; ++ci) {
          ierr = DMLabelSetValue(subpointmap, closure[2*ci], tempMark); CHKERRQ(ierr);
        }
        ierr = DMPlexRestoreTransitiveClosure(dm, point, PETSC_TRUE, &closureSize, &closure); CHKERRQ(ierr);
      }
    }
    /* Mark adjacent points */
    PetscBool updated = PETSC_TRUE;
    while (updated) {
      updated = PETSC_FALSE;
      for (p = 0; p < npoints; ++p) {
        point = points[p];
        ierr = DMLabelGetValue(subpointmap, point, &leafpoint); CHKERRQ(ierr);
        if (leafpoint == defaultValue) {
          /* Check if this point is adjacent to any already marked point */
          PetscInt val;
          const PetscInt *cone;
          ierr = DMPlexGetConeSize(dm, point, &closureSize); CHKERRQ(ierr);
          ierr = DMPlexGetCone(dm, point, &cone); CHKERRQ(ierr);
          for (ci = 0; ci < closureSize; ++ci) {
            ierr = DMLabelGetValue(subpointmap, cone[ci], &val); CHKERRQ(ierr);
            if (val != defaultValue) break;
          }
          if (ci == closureSize) continue;
          ierr = DMPlexGetTransitiveClosure(dm, point, PETSC_TRUE, &closureSize, &closure); CHKERRQ(ierr);
          for (ci = 0; ci < closureSize; ++ci) {
            ierr = DMLabelGetValue(subpointmap, closure[2*ci], &val); CHKERRQ(ierr);
            if (val == defaultValue) {
              ierr = DMLabelSetValue(subpointmap, closure[2*ci], tempMarkGhost); CHKERRQ(ierr);
            }
          }
          ierr = DMPlexRestoreTransitiveClosure(dm, point, PETSC_TRUE, &closureSize, &closure); CHKERRQ(ierr);
          updated = PETSC_TRUE;
        }
      }
    }
    ierr = ISRestoreIndices(marked, &points); CHKERRQ(ierr);
    ierr = ISDestroy(&marked); CHKERRQ(ierr);
  }
  PetscHMapIDestroy(&pointmap);

  /* Add overlap */
  {
    PetscInt overlap;
    if (addOverlap) {
        ierr = DMPlexGetOverlap(dm, &overlap); CHKERRQ(ierr);
        printf("overlap = %d\n", overlap);
    }
    else {
        overlap = 0;
    }
    while (overlap > 0) {
      ierr = DMLabelGetStratumIS(subpointmap, tempMark, &marked); CHKERRQ(ierr);
      if (marked) {
        PetscInt *adj = NULL;
        PetscInt  apoint, apointval; 
        ierr = ISGetLocalSize(marked, &npoints); CHKERRQ(ierr);
        ierr = ISGetIndices(marked, &points); CHKERRQ(ierr);
        for (p = 0; p < npoints; ++p) {
          PetscInt adjSize = PETSC_DETERMINE, a;
          point = points[p];
          ierr = DMPlexGetAdjacency(dm, point, &adjSize, &adj); CHKERRQ(ierr);
          for (a = 0; a < adjSize; ++a) {
            apoint = adj[a];
            ierr = DMLabelGetValue(subpointmap, apoint, &apointval); CHKERRQ(ierr);
            if (apointval == tempMarkGhost) {
              ierr = DMLabelSetValue(subpointmap, apoint, tempMark); CHKERRQ(ierr);
            }
          }
          PetscFree(adj);
        }
        ierr = ISRestoreIndices(marked, &points); CHKERRQ(ierr);
        ierr = ISDestroy(&marked); CHKERRQ(ierr);
      }
      overlap--;
    }
  }

  /* Remove redundant marks */
  ierr = DMLabelGetStratumIS(subpointmap, tempMarkGhost, &marked); CHKERRQ(ierr);
  if (marked) {
    ierr = ISGetLocalSize(marked, &npoints); CHKERRQ(ierr);
    ierr = ISGetIndices(marked, &points); CHKERRQ(ierr);
    for (p = 0; p < npoints; ++p) {
      ierr = DMLabelSetValue(subpointmap, points[p], defaultValue); CHKERRQ(ierr);
      }
    ierr = ISRestoreIndices(marked, &points); CHKERRQ(ierr);
    ierr = ISDestroy(&marked); CHKERRQ(ierr);
  }

  /* Mark subpointMap */
  ierr = PetscMalloc1(tdim+1, &pStart); CHKERRQ(ierr);
  ierr = PetscMalloc1(tdim+1, &pEnd); CHKERRQ(ierr);
  for (d = 0; d <= tdim; ++d) {
    ierr = DMPlexGetDepthStratum(dm, d, &pStart[d], &pEnd[d]);CHKERRQ(ierr);
  }
  ierr = DMLabelGetStratumIS(subpointmap, tempMark, &marked); CHKERRQ(ierr);
  if (marked) {
    ierr = ISGetLocalSize(marked, &npoints); CHKERRQ(ierr);
    ierr = ISGetIndices(marked, &points); CHKERRQ(ierr);
    for (p = 0; p < npoints; ++p) {
      point = points[p];
      for (d = 0; d <= tdim; ++d) {
        if ((point >= pStart[d]) && (point < pEnd[d])) {
          ierr = DMLabelSetValue(subpointmap, point, d); CHKERRQ(ierr);
          break;
        }
      }
    }
    ierr = ISRestoreIndices(marked, &points); CHKERRQ(ierr);
    ierr = ISDestroy(&marked); CHKERRQ(ierr);
  }
  ierr = PetscFree(pStart); CHKERRQ(ierr);
  ierr = PetscFree(pEnd); CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


static PetscErrorCode _DMPlexSubmeshSetPointSF(DM dm, DM subdm,
                                              const PetscInt *stratumOffsets,
                                              const PetscInt *stratumSizes,
                                              const PetscInt **stratumIndices)
{
  PetscErrorCode     ierr;
  PetscSF            sf, subsf;
  PetscMPIInt        rank, size, subsize;
  PetscInt           pStart, pEnd, subStart, subEnd, p, idx, point;//, subpoint;
  PetscInt           nroots, nleaves, nsubleaves, nsubpoints;
  IS                 subpointIS;
  const PetscInt    *subpoints;
  const PetscInt    *ilocal;
  const PetscSFNode *iremote;
  PetscSFNode       *subiremote;
  PetscInt          *subilocal;
  MPI_Comm           comm, subcomm;
  PetscInt           subtdim;
  PetscInt          *updatedOwnersReduced;
  PetscInt          *updatedOwners;
  PetscInt          *updatedIndicesReduced;
  PetscInt          *updatedIndices;

  PetscFunctionBegin;

  comm = PetscObjectComm((PetscObject)dm);
  subcomm = PetscObjectComm((PetscObject)subdm);
  ierr = MPI_Comm_rank(comm, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm, &size); CHKERRQ(ierr);
  ierr = MPI_Comm_size(subcomm, &subsize); CHKERRQ(ierr);
  if (subsize != size && subsize != 1) {
    SETERRQ(subcomm, PETSC_ERR_PLIB, "Not implemented for subcomm not of size 1");
  }
  ierr = DMGetPointSF(dm, &sf); CHKERRQ(ierr);
  ierr = DMGetPointSF(subdm, &subsf); CHKERRQ(ierr);
  if (subsize == 1) {
    PetscFunctionReturn(0);
  } 
  ierr = PetscSFGetGraph(sf, &nroots, &nleaves, &ilocal, &iremote); CHKERRQ(ierr);
  if (nroots < 0) {
    ierr = PetscSFSetGraph(subsf, 0, 0, NULL, PETSC_OWN_POINTER, NULL, PETSC_OWN_POINTER); CHKERRQ(ierr);
    PetscFunctionReturn(0);
  }
  ierr = DMPlexGetChart(dm, &pStart, &pEnd); CHKERRQ(ierr);
  if (nroots != (pEnd - pStart)) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Number of roots does not match chart");
  ierr = DMPlexGetChart(subdm, &subStart, &subEnd); CHKERRQ(ierr);
  ierr = DMGetDimension(subdm, &subtdim); CHKERRQ(ierr);

  ierr = DMPlexCreateSubpointIS(subdm, &subpointIS); CHKERRQ(ierr);
  if (subpointIS) {
    ierr = ISGetIndices(subpointIS, &subpoints); CHKERRQ(ierr);
    ierr = ISGetLocalSize(subpointIS, &nsubpoints); CHKERRQ(ierr);
  } else {
    nsubpoints = 0;
  }
  if (nsubpoints != (subEnd - subStart)) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Subpoint IS does not cover complete chart");

  /* update owners */
  ierr = PetscMalloc1(nroots, &updatedOwnersReduced); CHKERRQ(ierr);
  ierr = PetscMalloc1(nroots, &updatedOwners); CHKERRQ(ierr);
  for (p = pStart; p < pEnd; ++p) {
    updatedOwnersReduced[p-pStart] = -1;
    updatedOwners[p-pStart] = -1;
  }
  for (p = 0; p < nsubpoints; ++p) {
    updatedOwners[subpoints[p]] = rank;
  }

  ierr = PetscSFReduceBegin(sf, MPIU_INT, updatedOwners, updatedOwnersReduced, MPI_MAX); CHKERRQ(ierr);
  ierr = PetscSFReduceEnd(sf, MPIU_INT, updatedOwners, updatedOwnersReduced, MPI_MAX); CHKERRQ(ierr);

  for (p = 0; p < nsubpoints; ++p) {
    updatedOwnersReduced[subpoints[p]] = rank;
    updatedOwners[subpoints[p]] = rank;
  }

  ierr = PetscSFBcastBegin(sf, MPIU_INT, updatedOwnersReduced, updatedOwners); CHKERRQ(ierr);
  ierr = PetscSFBcastEnd(sf, MPIU_INT, updatedOwnersReduced, updatedOwners); CHKERRQ(ierr);

  ierr = PetscFree(updatedOwnersReduced); CHKERRQ(ierr);

  /* broadcast new remote subdm indices including ones on new owners */
  ierr = PetscMalloc1(nroots, &updatedIndicesReduced); CHKERRQ(ierr);
  ierr = PetscMalloc1(nroots, &updatedIndices); CHKERRQ(ierr);

  for (p = pStart; p < pEnd; ++p) {
    updatedIndicesReduced[p - pStart] = -1;
    updatedIndices[p - pStart] = -1;
  }
  for (p = 0; p < nsubpoints; ++p) {
    point = subpoints[p];
    if (updatedOwners[point] == rank) {
      updatedIndices[point] = p;
    }
  }
  ierr = PetscSFReduceBegin(sf, MPIU_INT, updatedIndices, updatedIndicesReduced, MPI_MAX); CHKERRQ(ierr);
  ierr = PetscSFReduceEnd(sf, MPIU_INT, updatedIndices, updatedIndicesReduced, MPI_MAX); CHKERRQ(ierr);

  for (p = 0; p < nsubpoints; ++p) {
    point = subpoints[p];
    updatedIndicesReduced[point] = p;
  }
  ierr = PetscSFBcastBegin(sf, MPIU_INT, updatedIndicesReduced, updatedIndices); CHKERRQ(ierr);
  ierr = PetscSFBcastEnd(sf, MPIU_INT, updatedIndicesReduced, updatedIndices); CHKERRQ(ierr);
  ierr = PetscFree(updatedIndicesReduced); CHKERRQ(ierr);

  /* set subsf */
  nsubleaves = 0;
  for (p = 0; p < nsubpoints; ++p) {
    if (updatedOwners[subpoints[p]] != rank) ++nsubleaves;
  }
  ierr = PetscMalloc1(nsubleaves, &subilocal); CHKERRQ(ierr);
  ierr = PetscMalloc1(nsubleaves, &subiremote); CHKERRQ(ierr);

  for (p = 0, idx = 0; p < nsubpoints; ++p) {
    if (updatedOwners[subpoints[p]] != rank) {
      subilocal[idx] = p;
      subiremote[idx].rank = updatedOwners[subpoints[p]];
      subiremote[idx].index = updatedIndices[subpoints[p]];
      if (subiremote[idx].rank < 0) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Invalid remote rank");
      if (subiremote[idx].index < 0) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Invalid remote subpoint");
      ++idx;
    }
  }
  if (subpointIS) {
    ierr = ISRestoreIndices(subpointIS, &subpoints); CHKERRQ(ierr);
    ierr = ISDestroy(&subpointIS); CHKERRQ(ierr);
  }
  //PetscHMapIDestroy(&pointmap);
  ierr = PetscFree(updatedOwners); CHKERRQ(ierr);
  ierr = PetscFree(updatedIndices); CHKERRQ(ierr);
  ierr = PetscSFSetGraph(subsf, subEnd - subStart, nsubleaves, subilocal, PETSC_OWN_POINTER, subiremote, PETSC_OWN_POINTER); CHKERRQ(ierr);

  PetscFunctionReturn(0);
};


static PetscErrorCode _DMPlexCreateSubDMPlex(DM dm, DM *subdm, DMLabel filter, PetscInt filterValue, PetscInt height)
{
  PetscErrorCode   ierr;
  MPI_Comm         comm, subcomm;
  IS               reordering;
  DMLabel          subpointMap;
  PetscInt         tdim, gdim;
  PetscInt         subtdim;
  PetscInt        *stratumSizes, *stratumOffsets;
  const PetscInt **stratumIndices;
  IS              *stratumISes;
  PetscInt         d;

  comm = PetscObjectComm((PetscObject)dm);
  subcomm = comm;
  ierr = DMCreate(subcomm, subdm);CHKERRQ(ierr);
  ierr = DMSetType(*subdm, DMPLEX);CHKERRQ(ierr);

  ierr = DMGetDimension(dm, &tdim);CHKERRQ(ierr);
  subtdim = tdim - height;
  if (subtdim < 0) {
    SETERRQ(comm, PETSC_ERR_PLIB, "Cannot create mesh with negative topological dimension");
  }
  ierr = DMSetDimension(*subdm, subtdim);CHKERRQ(ierr);
  ierr = DMGetCoordinateDim(dm, &gdim);CHKERRQ(ierr);
  ierr = DMSetCoordinateDim(*subdm, gdim);CHKERRQ(ierr);

  /* Copy information relevant to adjacency */
  {
    PetscBool   useCone, useClosure;
    ierr = DMGetAdjacency(dm, -1, &useCone, &useClosure); CHKERRQ(ierr);
    ierr = DMSetAdjacency(*subdm, -1, useCone, useClosure); CHKERRQ(ierr);
  }
  ierr = DMCopyFields(dm, *subdm);CHKERRQ(ierr);
  {
    PetscBool       useAnchors;
    PetscErrorCode (*useradjacency)(DM,PetscInt,PetscInt*,PetscInt[],void*);
    void           *useradjacencyctx;
    ierr = DMPlexGetAdjacencyUser(dm, &useradjacency, &useradjacencyctx); CHKERRQ(ierr);
    ierr = DMPlexSetAdjacencyUser(*subdm, useradjacency, useradjacencyctx); CHKERRQ(ierr);
    ierr = DMPlexGetAdjacencyUseAnchors(dm, &useAnchors); CHKERRQ(ierr);
    ierr = DMPlexSetAdjacencyUseAnchors(*subdm, useAnchors); CHKERRQ(ierr);
  }

  /* Create non-overlapping subpointMap */
  {
    DMLabel            nonOverlappingFilter;
    PetscSF            sf;
    PetscInt           nroots, nleaves;
    const PetscInt    *ilocal;
    const PetscSFNode *iremote;
    IS                 marked;
    PetscInt           p, point, npoints, leafpoint, adjpoint;
    const PetscInt    *points;

    ierr = DMLabelCreate(PetscObjectComm((PetscObject)dm), "nonOverlappingFilter", &nonOverlappingFilter); CHKERRQ(ierr);
    ierr = DMGetPointSF(dm, &sf); CHKERRQ(ierr);
    ierr = PetscSFGetGraph(sf, &nroots, &nleaves, &ilocal, &iremote); CHKERRQ(ierr);
    ierr = DMLabelCreate(comm, "subpoint_map", &subpointMap); CHKERRQ(ierr);
    ierr = _DMPlexMarkSubpointMap_Closure(dm, filter, filterValue, height, subpointMap, PETSC_TRUE); CHKERRQ(ierr);
    ierr = DMLabelDestroy(&nonOverlappingFilter); CHKERRQ(ierr);
  }
  ierr = PetscMalloc1(subtdim+1, &stratumSizes); CHKERRQ(ierr);
  ierr = PetscMalloc1(subtdim+1, &stratumOffsets); CHKERRQ(ierr);
  ierr = PetscMalloc1(subtdim+1, &stratumIndices); CHKERRQ(ierr);
  ierr = PetscMalloc1(subtdim+1, &stratumISes); CHKERRQ(ierr);
  for (d = 0; d <= subtdim; ++d) {
    ierr = DMLabelGetStratumSize(subpointMap, d, &stratumSizes[d]); CHKERRQ(ierr);
    ierr = DMLabelGetStratumIS(subpointMap, d, &stratumISes[d]); CHKERRQ(ierr);
    if (stratumISes[d]) { 
        ierr = ISGetIndices(stratumISes[d], &stratumIndices[d]); CHKERRQ(ierr);
    }
  }

  stratumOffsets[subtdim] = 0;
  if (subtdim > 0) {
    stratumOffsets[0] = stratumOffsets[subtdim] + stratumSizes[subtdim];
  }
  if (subtdim > 1) {
    stratumOffsets[subtdim - 1] = stratumOffsets[0] + stratumSizes[0];
  }
  if (subtdim > 2) {
    stratumOffsets[subtdim - 2] = stratumOffsets[subtdim - 1] + stratumSizes[subtdim - 1];
  }
  if (subtdim > 3) {
    SETERRQ(subcomm, PETSC_ERR_PLIB, "Only coded for max 3 dimensional DMs");
  }

  ierr = DMPlexSetSubpointMap(*subdm, subpointMap);CHKERRQ(ierr);
  ierr = DMLabelDestroy(&subpointMap); CHKERRQ(ierr);

  ierr = DMPlexSubmeshSetTopology(dm, *subdm, stratumOffsets, stratumSizes, stratumIndices); CHKERRQ(ierr);
  ierr = DMPlexSubmeshSetCoordinates(dm, *subdm, stratumOffsets, stratumSizes, stratumIndices); CHKERRQ(ierr);
  ierr = _DMPlexSubmeshSetPointSF(dm, *subdm, stratumOffsets, stratumSizes, stratumIndices); CHKERRQ(ierr);

  /* Add the partition overlap to the subdm */
  /*{
    DM              dmOverlap;
    PetscSF         sfOverlap;
    const PetscInt  overlap = 1;//((DM_Plex*)dm->data)->overlap;
    ierr = DMPlexDistributeOverlap(*subdm, overlap, &sfOverlap, &dmOverlap); CHKERRQ(ierr);
    ierr = DMDestroy(subdm); CHKERRQ(ierr);
    *subdm = dmOverlap;
  }*/

  {
  PetscSF            subsf;
  PetscInt           nroots, nleaves;
  //IS                 subpointIS;
  //const PetscInt    *subpoints;
  const PetscInt    *ilocal;
  const PetscSFNode *iremote;
  PetscInt           rank, size;//, nsubpoints;
  //const PetscInt     p;

  ierr = MPI_Comm_rank(comm, &rank); CHKERRQ(ierr);
  ierr = MPI_Comm_size(comm, &size); CHKERRQ(ierr);

  ierr = DMGetPointSF(*subdm, &subsf); CHKERRQ(ierr);
  ierr = PetscSFGetGraph(subsf, &nroots, &nleaves, &ilocal, &iremote); CHKERRQ(ierr);
  }



  /* Finalize */
  for (d = 0; d <= subtdim; d++) {
    if (stratumISes[d]) { 
        ierr = ISRestoreIndices(stratumISes[d], &stratumIndices[d]); CHKERRQ(ierr);
    }
    ierr = ISDestroy(&stratumISes[d]); CHKERRQ(ierr);
  }
  ierr = PetscFree(stratumSizes); CHKERRQ(ierr);
  ierr = PetscFree(stratumOffsets); CHKERRQ(ierr);
  ierr = PetscFree(stratumIndices); CHKERRQ(ierr);
  ierr = PetscFree(stratumISes); CHKERRQ(ierr);

  ierr = DMPlexCheckFaces(*subdm, 0); CHKERRQ(ierr);
  ierr = DMPlexCheckSymmetry(*subdm); CHKERRQ(ierr);
  ierr = DMPlexGetOrdering(*subdm, MATORDERINGRCM, NULL, &reordering); CHKERRQ(ierr);

  ierr = ISDestroy(&reordering); CHKERRQ(ierr);

  PetscFunctionReturn(0);
};






/*

test 0:
Extract left half domain (2 proc.)
Check ownership transferring

Global numbering:

           6--16---7---17--11--22--12--23--13
           |       |       |       |       |
          18   0   19   1  24  2   25   3  26
           |       |       |       |       |
           4--14---5---15--8---20--9---21--10

Local numberings:

                            dm                                        subdm

           5--10---6---11-(7)                           5--10---6---11--7
           |       |       |                            |       |       |
rank 0:   12   0   13  1  (14)                    -->  12   0   13  1   14
           |       |       |                            |       |       |
           2---8---3---9--(4)                           2---8---3---9---4


                           5--10---6---11---7
                           |       |        |
rank 1:                   12   0   13   1   14    -->   None
                           |       |        |
                           2---8---3----9---4



-overlap 1:

sub cells = {0, 1, 3} (global)

* ownership transfer
* disconnectedness


           5--13---6--14--(9)-(18)(10)                  4--10---5---11--7
           |       |       |       |                    |       |       |
rank 0:   15   0   16  1  (19) (2)(20)                 12   0   13  1   14
           |       |       |       |                    |       |       |
           3--11---4--12--(7)-(17)(8)                   2---8---3---9---6

                 (10)-(19)-6--13---7---14--8                                    3---6---4
                   |       |       |       |                                    |       |
rank 1:          (20) (2)  15  0   16  1   17                                   7   0   8
                   |       |       |       |                                    |       |
                  (9)-(18)-3--11---4---12--5                                    1 --5---2




*/

typedef struct {
  PetscInt  testNum;                      /* Indicates the mesh to create */
  PetscInt  overlap;                      /* The partition overlap */
} AppCtx;


PetscErrorCode ProcessOptions(MPI_Comm comm, AppCtx *options)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  options->testNum      = 0;
  options->overlap      = 0;

  ierr = PetscOptionsBegin(comm, "", "Meshing Interpolation Test Options", "DMPLEX");CHKERRQ(ierr);
  ierr = PetscOptionsBoundedInt("-testNum", "The mesh to create", "ex37.c", options->testNum, &options->testNum, NULL,0);CHKERRQ(ierr);
  ierr = PetscOptionsBoundedInt("-overlap", "The partition overlap", "ex37.c", options->overlap, &options->overlap, NULL,0);CHKERRQ(ierr);
  ierr = PetscOptionsEnd();
  PetscFunctionReturn(0);
}


int main(int argc, char **argv)
{
  DM               dm, subdm;
  DMLabel          filter;
  PetscInt         height;
  const PetscInt   filterValue = 1;
  MPI_Comm         comm;
  PetscMPIInt rank, size;
  AppCtx           user;
  PetscErrorCode   ierr;

  ierr = PetscInitialize(&argc, &argv, NULL, NULL); if (ierr) return ierr;
  comm = PETSC_COMM_WORLD;
  ierr = ProcessOptions(comm, &user);CHKERRQ(ierr);

  ierr = MPI_Comm_rank(comm, &rank); CHKERRQ(ierr);
  ierr = DMLabelCreate(comm, "filter", &filter); CHKERRQ(ierr);

  switch (user.testNum) {
    case 0:
    {
      /* Create parallel dm */
      const PetscInt faces[2] = {4,1};
      DM             pdm;
      PetscSF        sf;
      ierr = DMPlexCreateBoxMesh(comm, 2, PETSC_FALSE, faces, NULL, NULL, NULL, PETSC_TRUE, &dm); CHKERRQ(ierr);
      ierr = DMPlexDistribute(dm, user.overlap, &sf, &pdm); CHKERRQ(ierr);
      if (pdm) {
        ierr = DMDestroy(&dm); CHKERRQ(ierr);
        dm = pdm;
      }
      if (sf) {
        ierr = PetscSFDestroy(&sf); CHKERRQ(ierr);
      }
      /* Create filter label */
      height = 0;
      switch (user.overlap) {
        case 0:
          if (rank==0) {
            DMLabelSetValue(filter, 0, filterValue);
            DMLabelSetValue(filter, 1, filterValue);
          }
          break;
        case 1:
          if (rank==0) {
            DMLabelSetValue(filter, 0, filterValue);
            DMLabelSetValue(filter, 1, filterValue);
          }
          if (rank==1) {
            DMLabelSetValue(filter, 2, filterValue);
            DMLabelSetValue(filter, 1, filterValue);
          }
          break;
        case 2:
          break;
      }
      break;
    }
  }
// testing........
{
  PetscInt pStart, pEnd, p;
  ierr = DMPlexGetChart(dm, &pStart, &pEnd); CHKERRQ(ierr);
  IS        globalPointNumbers;
  const PetscInt *ixs;
  ierr = DMPlexCreatePointNumbering(dm, &globalPointNumbers);CHKERRQ(ierr);
  ierr = ISGetIndices(globalPointNumbers, &ixs);CHKERRQ(ierr);
  for (p = pStart; p < pEnd; ++p) {
    printf("rank = %d, (globalnum, localnum) = (%d, %d)\n", rank, ixs[p-pStart], p);
  }
  ierr = ISRestoreIndices(globalPointNumbers, &ixs);CHKERRQ(ierr);
  ierr = ISDestroy(&globalPointNumbers);CHKERRQ(ierr);
}
// testing end

  ierr = _DMPlexCreateSubDMPlex(dm, &subdm, filter, filterValue, height);
  ierr = DMLabelDestroy(&filter); CHKERRQ(ierr);

  ierr = PetscObjectSetName((PetscObject) dm, "Example_DM");CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) subdm, "Example_SubDM");CHKERRQ(ierr);
  ierr = DMViewFromOptions(dm, NULL, "-dm_view"); CHKERRQ(ierr);
  ierr = DMViewFromOptions(subdm, NULL, "-dm_view"); CHKERRQ(ierr);

  ierr = DMDestroy(&dm); CHKERRQ(ierr);
  ierr = DMDestroy(&subdm); CHKERRQ(ierr);
  ierr = PetscFinalize();
  return ierr;
}

/*TEST

  # Two cell test meshes 0-7
  test:
    suffix: 0
    nsize: 2
    args: -testNum 0 -overlap 0 -dm_view ascii::ascii_info_detail
  test:
    suffix: 1
    nsize: 2
    args: -testNum 0 -overlap 1 -dm_view ascii::ascii_info_detail
  test:
    suffix: 2
    args: -dim 2 -cell_simplex 0 -dm_view ascii::ascii_info_detail
  test:
    suffix: 3
    nsize: 2
    args: -petscpartitioner_type simple -dim 2 -cell_simplex 0 -dm_view ascii::ascii_info_detail

TEST*/
