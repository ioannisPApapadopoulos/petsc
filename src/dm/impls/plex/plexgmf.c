#include <petsc/private/dmpleximpl.h>   /*I      "petscdmplex.h"   I*/
#include "libmesh7.h"



PetscErrorCode DMPlexWrite_gmfMesh2d(DM dm, const char bdLabelName[], const char meshName[], const PetscBool ascii) {
  
  DM                 cdm;
  PetscSection       coordSection;
  Vec                coordinates;
  const PetscScalar *coords;
  PetscInt           dim, cStart, cEnd, numCells, c, vStart, vEnd, numVertices, v, off;
  PetscInt           idx[3], i;
  PetscBool          B64;
  char               fileName[512];
  long long          meshIndex;
  int                fileVersion;
  PetscErrorCode     ierr;
  
  
  
  B64 = PETSC_TRUE;
  ierr = DMGetDimension(dm, &dim);CHKERRQ(ierr);
  ierr = DMGetDimension(dm, &dim);CHKERRQ(ierr);
  ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd);CHKERRQ(ierr);
  numCells = cEnd - cStart;
  ierr = DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd);CHKERRQ(ierr);
  numVertices = vEnd - vStart;
  
  strcpy(fileName, meshName);

  if ( ascii )
    strcat(fileName, ".mesh");
  else
    strcat(fileName, ".meshb");
  
  if ( B64 )
    fileVersion = GmfDouble;
  else
    fileVersion = GmfFloat; 
     

  //--- Open file
  if ( !(meshIndex = GmfOpenMesh(fileName,GmfWrite,fileVersion,dim)) ) {
    fprintf(stderr,"####  ERROR: mesh file %s cannot be opened\n", fileName);
    exit(1);
  }
  printf("  %%%% %s opened\n",fileName);
  
  
  
  GmfSetKwd(meshIndex, GmfVertices, numVertices);
  ierr = DMGetCoordinateDM(dm, &cdm);CHKERRQ(ierr);
  ierr = DMGetDefaultSection(cdm, &coordSection);CHKERRQ(ierr);
  ierr = DMGetCoordinatesLocal(dm, &coordinates);CHKERRQ(ierr);
  ierr = VecGetArrayRead(coordinates, &coords);CHKERRQ(ierr);
  for (v = vStart; v < vEnd; ++v) {
    ierr = PetscSectionGetOffset(coordSection, v, &off);CHKERRQ(ierr);
    GmfSetLin(meshIndex, GmfVertices, coords[off], coords[off+1], 0);  
  }
  ierr = VecRestoreArrayRead(coordinates, &coords);CHKERRQ(ierr);
  
  GmfSetKwd(meshIndex, GmfTriangles, numCells);
  for (c = cStart; c < cEnd; ++c) {
//    printf("DEBUG   Triangle %d\n", c-cStart);
    PetscInt *closure = NULL;
    PetscInt  closureSize, cl, p;
    ierr = DMPlexGetTransitiveClosure(dm, c, PETSC_TRUE, &closureSize, &closure);CHKERRQ(ierr);
    for (cl = 0, i=0; cl < closureSize*2; cl += 2) {
      p = closure[cl];
//      printf("DEBUG    cl: %d\n", cl);
      if (p >= vStart && p < vEnd) {
//        printf("DEBUG     idx[%d] = %d\n", c-cStart, i, cl-vStart+1);
        idx[i++] = p-vStart+1;
      }
    }
    GmfSetLin(meshIndex, GmfTriangles, idx[0], idx[1], idx[2], 0);  
    ierr = DMPlexRestoreTransitiveClosure(dm, c, PETSC_TRUE, &closureSize, &closure);CHKERRQ(ierr);
  }
  
  if ( !GmfCloseMesh(meshIndex) ) {
    fprintf(stderr,"#### ERROR: mesh file %s cannot be closed\n",fileName);
    exit(1);
  }

  PetscFunctionReturn(0);
  
}