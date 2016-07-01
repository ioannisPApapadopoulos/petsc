#include <petsc/private/dmpleximpl.h>   /*I      "petscdmplex.h"   I*/
#include "libmesh7.h"


// solTypes: 1 for a scalar, 2 for a vector, 3 for symmetric matrix and 4 for a full matrix

PetscErrorCode DMPlexWrite_gmfMesh2d(DM dm, PetscBool writeMesh, PetscInt numSol, Vec * sol,  PetscInt * solTypes, 
                                      const char bdLabelName[], const char meshName[], const char * solNames[], 
                                      const PetscBool ascii) {
  
  DM                  cdm;
  PetscSection        coordSection;
  Vec                 coordinates;
  const PetscScalar * coords, * solution;
  PetscScalar       * buffer;
  PetscInt            dim, cStart, cEnd, numCells, c, vStart, vEnd, numVertices, v, off;
  PetscInt            idx[3], i, iSol;
  PetscBool           B64;
  char                fileName[512];
  long long           meshIndex, solIndex;
  int                 fileVersion, solKeyword;
  PetscErrorCode      ierr;
  
  
  
  B64 = PETSC_TRUE;
  ierr = DMGetDimension(dm, &dim);CHKERRQ(ierr);
  ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd);CHKERRQ(ierr);
  numCells = cEnd - cStart;
  ierr = DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd);CHKERRQ(ierr);
  numVertices = vEnd - vStart;
  ierr = DMGetCoordinateDM(dm, &cdm);CHKERRQ(ierr);
  ierr = DMGetDefaultSection(cdm, &coordSection);CHKERRQ(ierr);


  if (writeMesh){

    strcpy(fileName, meshName);
    if ( ascii )
      strcat(fileName, ".mesh");
    else
      strcat(fileName, ".meshb");    
    if ( B64 )
      fileVersion = GmfDouble;
    else
      fileVersion = GmfFloat; 
       
    if ( !(meshIndex = GmfOpenMesh(fileName, GmfWrite, fileVersion, dim)) ) {
      fprintf(stderr,"####  ERROR: mesh file %s cannot be opened\n", fileName);
      exit(1);
    }
    printf("  %%%% %s opened\n",fileName);
    
    GmfSetKwd(meshIndex, GmfVertices, numVertices);
    ierr = DMGetCoordinatesLocal(dm, &coordinates);CHKERRQ(ierr);
    ierr = VecGetArrayRead(coordinates, &coords);CHKERRQ(ierr);
    for (v = vStart; v < vEnd; ++v) {
      ierr = PetscSectionGetOffset(coordSection, v, &off);CHKERRQ(ierr);
      GmfSetLin(meshIndex, GmfVertices, coords[off], coords[off+1], 0);  
    }
    ierr = VecRestoreArrayRead(coordinates, &coords);CHKERRQ(ierr);
    
    GmfSetKwd(meshIndex, GmfTriangles, numCells);
    for (c = cStart; c < cEnd; ++c) {
      PetscInt *closure = NULL;
      PetscInt  closureSize, cl, p;
      ierr = DMPlexGetTransitiveClosure(dm, c, PETSC_TRUE, &closureSize, &closure);CHKERRQ(ierr);
      for (cl = 0, i=0; cl < closureSize*2; cl += 2) {
        p = closure[cl];
        if (p >= vStart && p < vEnd) 
          idx[i++] = p-vStart+1;
      }
      GmfSetLin(meshIndex, GmfTriangles, idx[0], idx[1], idx[2], 0);  
      ierr = DMPlexRestoreTransitiveClosure(dm, c, PETSC_TRUE, &closureSize, &closure);CHKERRQ(ierr);
    }
    
    if ( !GmfCloseMesh(meshIndex) ) {
      fprintf(stderr,"#### ERROR: mesh file %s cannot be closed\n",fileName);
      exit(1);
    }
  }


  for (iSol = 0; iSol < numSol; ++iSol) {

    // TODO: a terme il faudra gerer les sol non scalaire et pouvoir choisir entre toutes dans un meme fichier ou pas

    solKeyword = GmfSolAtVertices;
    strcpy(fileName, solNames[iSol]);
    if ( ascii )
      strcat(fileName, ".sol");
    else
      strcat(fileName, ".solb");    
    if ( B64 )
      fileVersion = GmfDouble;
    else
      fileVersion = GmfFloat; 
       
    if ( !(solIndex = GmfOpenMesh(fileName, GmfWrite, fileVersion, dim)) ) {
      fprintf(stderr,"####  ERROR: solution file %s cannot be opened\n", fileName);
      exit(1);
    }
    printf("  %%%% %s opened\n",fileName);

    GmfSetKwd(solIndex, solKeyword, numVertices, 1, &solTypes[iSol]);
    ierr = VecGetArrayRead(sol[iSol], &solution);CHKERRQ(ierr);
    switch(solTypes[iSol]) {
      case 1 :
        ierr = PetscMalloc1(1, &buffer);CHKERRQ(ierr);
        break;
      case 2 :
        ierr = PetscMalloc1(2, &buffer);CHKERRQ(ierr);
        break;
      default :
        printf("####  ERROR  non-scalar solutions not implemented yet\n");
        exit(1);
    }
    for (v = vStart; v < vEnd; ++v) {
      ierr = PetscSectionGetOffset(coordSection, v, &off);CHKERRQ(ierr);
      switch(solTypes[iSol]) {
        case 1 :
          buffer[0] = solution[off/dim];
          break;
        case 2 :
          buffer[0] = solution[off];
          buffer[1] = solution[off+1];
          break;
        default :
          printf("####  ERROR  non-scalar solutions not implemented yet\n");
          exit(1);
      }
      GmfSetLin(solIndex, solKeyword, buffer);      
    }
    ierr = VecRestoreArrayRead(sol[iSol], &solution);CHKERRQ(ierr);

    ierr = PetscFree(buffer);CHKERRQ(ierr);
    if ( !GmfCloseMesh(meshIndex) ) {
      fprintf(stderr,"#### ERROR: mesh file %s cannot be closed\n",fileName);
      exit(1);
    }

  }
  


  PetscFunctionReturn(0);
  
}



