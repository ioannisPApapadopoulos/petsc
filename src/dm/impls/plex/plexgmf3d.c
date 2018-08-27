#include <petsc/private/dmpleximpl.h>   /*I      "petscdmplex.h"   I*/
#include "libmesh6.h"





#undef __FUNCT__
#define __FUNCT__ "DMPlexWrite_gmfMesh3d_1sol"
// Temporary - For simple python use
PetscErrorCode DMPlexWrite_gmfMesh3d_1sol(DM dm, PetscBool writeMesh, const char bdLabelName[], const char meshName[], 
                                          Vec sol,  PetscInt solType, const char solNames[], 
                                          PetscSection section, PetscBool ascii) {
  
  PetscErrorCode ierr;
        
  PetscFunctionBegin;
  ierr = DMPlexWrite_gmfMesh3d(dm, writeMesh, bdLabelName, meshName, 1, &sol,  &solType, &solNames, section, ascii);CHKERRQ(ierr);                                
  PetscFunctionReturn(0);                                      
}


#undef __FUNCT__
#define __FUNCT__ "DMPlexWrite_gmfMesh3d_noSol"
// Temporary - For simple python use
PetscErrorCode DMPlexWrite_gmfMesh3d_noSol(DM dm, const char bdLabelName[], const char meshName[], PetscSection section, PetscBool ascii) {
  
  PetscErrorCode ierr;
        
  PetscFunctionBegin;
  ierr = DMPlexWrite_gmfMesh3d(dm, PETSC_TRUE, bdLabelName, meshName, 0, NULL, NULL, NULL, section, ascii);CHKERRQ(ierr);                                
  PetscFunctionReturn(0);                                      
}




#undef __FUNCT__
#define __FUNCT__ "DMPlexWrite_gmfMesh3d"
// solTypes: 1 for a scalar, 2 for a vector, 3 for symmetric matrix (only 6 values given: upper triangular matrix) and 4 for a full matrix
//            I add solType = 5 for symmetric matrixes given as 3x3 matrix

PetscErrorCode DMPlexWrite_gmfMesh3d(DM dm, PetscBool writeMesh, const char bdLabelName[], const char meshName[], 
                                     PetscInt numSol, Vec * sol,  PetscInt * solTypes, const char * solNames[],
                                     PetscSection section, PetscBool ascii) {
  
  DM                  cdm;
  PetscSection        coordSection;
  Vec                 coordinates;
  DMLabel             bdLabel = NULL;
  IS                  bdLabelIds;
  const PetscScalar   * coords, * solution;
  PetscScalar         * buffer;
  PetscInt            dim, cStart, cEnd, numCells, c, vStart, vEnd, numVertices, v, fStart, fEnd, numFacets, f, off;
  PetscInt            idx[4], i, iSol, tag, bdLabelSize, size, closureSize, cl, p;
  PetscInt            * closure;
  const PetscInt      * bdLabelVal;
  PetscBool           B64, flg=PETSC_FALSE;
  char                fileName[512];
  long long           meshIndex, solIndex;
  int                 fileVersion, solKeyword;
  PetscErrorCode      ierr;
  
  
  PetscFunctionBegin;
  B64 = PETSC_TRUE;
  ierr = DMGetDimension(dm, &dim);CHKERRQ(ierr);
  ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd);CHKERRQ(ierr);
  numCells = cEnd - cStart;
  ierr = DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd);CHKERRQ(ierr);
  numVertices = vEnd - vStart;
  ierr = DMPlexGetHeightStratum(dm, 1, &fStart, &fEnd);CHKERRQ(ierr);
  if (section) {
    coordSection = section;
  }
  else {
    ierr = DMGetCoordinateDM(dm, &cdm);CHKERRQ(ierr);
    ierr = DMGetDefaultSection(cdm, &coordSection);CHKERRQ(ierr);
  }

  if (writeMesh){
    strcpy(fileName, meshName);
    if ( ascii ) strcat(fileName, ".mesh");
    else         strcat(fileName, ".meshb");
    if ( B64 ) fileVersion = GmfDouble;
    else       fileVersion = GmfFloat; 
    if ( !(meshIndex = GmfOpenMesh(fileName, GmfWrite, fileVersion, dim)) ) {
      fprintf(stderr,"####  ERROR: mesh file %s cannot be opened\n", fileName);
      exit(1);
    }
    printf("  %%%% %s opened\n",fileName);
    GmfSetKwd(meshIndex, GmfVertices, numVertices);
    tag = 0;
    ierr = DMGetCoordinatesLocal(dm, &coordinates);CHKERRQ(ierr);
    ierr = VecGetArrayRead(coordinates, &coords);CHKERRQ(ierr);
    for (v = vStart; v < vEnd; ++v) {
      ierr = PetscSectionGetOffset(coordSection, v, &off);CHKERRQ(ierr);
      GmfSetLin(meshIndex, GmfVertices, coords[off], coords[off+1], coords[off+2], tag);  
    }
    ierr = VecRestoreArrayRead(coordinates, &coords);CHKERRQ(ierr);

    GmfSetKwd(meshIndex, GmfTetrahedra, numCells);
    tag = 0;
    for (c = cStart; c < cEnd; ++c) {
      closure = NULL;
      ierr = DMPlexGetTransitiveClosure(dm, c, PETSC_TRUE, &closureSize, &closure);CHKERRQ(ierr);
      for (cl = 0, i=0; cl < closureSize*2; cl += 2) {
        p = closure[cl];
        if (p >= vStart && p < vEnd) 
          idx[i++] = p-vStart+1;
      }
      GmfSetLin(meshIndex, GmfTetrahedra, idx[0], idx[1], idx[2], idx[3], tag);  
      ierr = DMPlexRestoreTransitiveClosure(dm, c, PETSC_TRUE, &closureSize, &closure);CHKERRQ(ierr);
    }


    if (bdLabelName) {
      ierr = PetscStrcmp(bdLabelName, "", &flg);CHKERRQ(ierr);
      if (!flg) {
        ierr = DMGetLabel(dm, bdLabelName, &bdLabel);CHKERRQ(ierr);
      }
    }
    if (bdLabel)  {
      ierr = DMGetLabelSize(dm, bdLabelName, &bdLabelSize);CHKERRQ(ierr);
      ierr = DMGetLabelIdIS(dm, bdLabelName, &bdLabelIds);CHKERRQ(ierr);
      ierr = ISGetIndices(bdLabelIds, &bdLabelVal);CHKERRQ(ierr);
      numFacets = 0;
      for (i=0; i<bdLabelSize; ++i){
        ierr = DMLabelGetStratumSize(bdLabel, bdLabelVal[i], &size);CHKERRQ(ierr);
        numFacets += size;
      }
      GmfSetKwd(meshIndex, GmfTriangles, numFacets);
      for (f = fStart; f < fEnd; ++f) {
      	ierr = DMLabelGetValue(bdLabel, f, &tag);CHKERRQ(ierr);
      	if (tag < 0) continue;
        closure = NULL;
        ierr = DMPlexGetTransitiveClosure(dm, f, PETSC_TRUE, &closureSize, &closure);CHKERRQ(ierr);
        if (closureSize != 7) {printf("####  ERROR  closure of a Triangle != 7, facet: %d, size: %d\n", f, closureSize);exit(1);}
        for (cl=0, i=0; cl<2*closureSize && i<3; cl+=2) {
          if (closure[cl] >= vStart && closure[cl] < vEnd ){
            idx[i] = closure[cl] - vStart + 1;
            i++;
          }
        }
        if (i != 3) {printf("####  ERROR  a Triangle has fewer than 3 vertices, facet: %d, size: %d\n", f, i);exit(1);}
        GmfSetLin(meshIndex, GmfTriangles, idx[0], idx[1], idx[2], tag);
        ierr = DMPlexRestoreTransitiveClosure(dm, f, PETSC_TRUE, &closureSize, &closure);CHKERRQ(ierr);
      }
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
    if ( ascii ) strcat(fileName, ".sol");
    else         strcat(fileName, ".solb");    
    if ( B64 ) fileVersion = GmfDouble;
    else       fileVersion = GmfFloat; 
       
    if ( !(solIndex = GmfOpenMesh(fileName, GmfWrite, fileVersion, dim)) ) {
      fprintf(stderr,"####  ERROR: solution file %s cannot be opened\n", fileName);
      exit(1);
    }
    printf("  %%%% %s opened\n",fileName);

    if (solTypes[iSol] == 5) { 
      solTypes[iSol] = 3;
      GmfSetKwd(solIndex, solKeyword, numVertices, 1, &solTypes[iSol]);
      solTypes[iSol] = 5;
    }
    else GmfSetKwd(solIndex, solKeyword, numVertices, 1, &solTypes[iSol]);
    
    ierr = VecGetArrayRead(sol[iSol], &solution);CHKERRQ(ierr);
    buffer = NULL;
    switch(solTypes[iSol]) {
      case 1 :
        ierr = PetscMalloc1(1, &buffer);CHKERRQ(ierr);
        break;
      case 2 :
        ierr = PetscMalloc1(3, &buffer);CHKERRQ(ierr);
        break;
      case 3 :
      case 5 :
        ierr = PetscMalloc1(6, &buffer);CHKERRQ(ierr);
        break;
      default :
        printf("####  ERROR  non-scalar solutions not implemented yet: %d\n", solTypes[iSol]);
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
          buffer[2] = solution[off+2];
          break;
        case 3 :
          off /= dim ;
          for (i=0; i<6; ++i)
          buffer[i] = solution[6*off+i];
          break;
        case 5 :
          off /= dim ;
          // be careful, libmesh reads the matrix by columns
          buffer[0] = solution[9*off]; buffer[1] = solution[9*off+1]; buffer[3] = solution[9*off+2];
                                       buffer[2] = solution[9*off+4]; buffer[4] = solution[9*off+5];
                                                                      buffer[5] = solution[9*off+8];
          break;
        default :
          printf("####  ERROR  non-scalar solutions not implemented yet!\n");
          exit(1);
      }
      GmfSetLin(solIndex, solKeyword, buffer);
    }
    ierr = VecRestoreArrayRead(sol[iSol], &solution);CHKERRQ(ierr);

    ierr = PetscFree(buffer);CHKERRQ(ierr);
    if ( !GmfCloseMesh(solIndex) ) {
      fprintf(stderr,"#### ERROR: solution file %s cannot be closed\n",fileName);
      exit(1);
    }

  }

  PetscFunctionReturn(0);
}





#undef __FUNCT__
#define __FUNCT__ "DMPlexCreateGmfFromFile_3d"
PetscErrorCode DMPlexCreateGmfFromFile_3d(const char meshName[], const char bdLabelName[], DM *dm)
{

  MPI_Comm            comm;
  PetscSection        coordSection;
  Vec                 coordinates;
  DMLabel             bdLabel = NULL;
  PetscScalar         * coordsIn, * coords;
  PetscScalar         buffer[3];
  PetscInt            bufferTet[4], bufferFac[3];
  PetscInt            dim, numCells, c, numVertices, v, numFacets, f, coordSize;
  PetscInt            i, tag, joinSize;
  const PetscInt      * join;
  PetscBool           flg;
  char                fileName[512];
  long long           meshIndex;
  int                 gmfVersion;
  PetscErrorCode      ierr;

  PetscFunctionBegin;
  comm = PETSC_COMM_WORLD;
  ierr = DMCreate(comm, dm);CHKERRQ(ierr);
  ierr = DMSetType(*dm, DMPLEX);CHKERRQ(ierr);
 
  strcpy(fileName, meshName);
  strcat(fileName, ".meshb");
  if ( !(meshIndex = GmfOpenMesh(fileName, GmfRead, &gmfVersion, &dim)) ) {
    strcpy(fileName, meshName);
    strcat(fileName,".mesh");
    if ( !(meshIndex = GmfOpenMesh(fileName, GmfRead, &gmfVersion, &dim)) ) {
      fprintf(stderr,"####  ERROR Mesh file %s.mesh[b] not found ", meshName);
      exit(1);
    }    
  }
  if (dim != 3) {
    printf("####  ERROR  Wrong dimension: %d != 3\n", dim);
    exit(1);
  }

  numVertices = GmfStatKwd(meshIndex, GmfVertices);
  numCells    = GmfStatKwd(meshIndex, GmfTetrahedra);

  if (numVertices <= 0 ) {
    printf("####  ERROR  Number of vertices: %d <= 0\n", numVertices);
    exit(1);
  }

  ierr = PetscMalloc1(numVertices*dim, &coordsIn);CHKERRQ(ierr);
  GmfGotoKwd(meshIndex, GmfVertices);
  for (v = 0; v < numVertices ; ++v) {
    GmfGetLin(meshIndex, GmfVertices, &buffer[0], &buffer[1], &buffer[2], &tag);
    coordsIn[dim*v]    = (double)buffer[0];
    coordsIn[dim*v+1]  = (double)buffer[1];
    coordsIn[dim*v+2]  = (double)buffer[2];
  }

  GmfGotoKwd(meshIndex, GmfTetrahedra);
  ierr = DMPlexSetChart(*dm, 0, numCells+numVertices);CHKERRQ(ierr);
  for (c = 0; c < numCells; ++c) {
    ierr = DMPlexSetConeSize(*dm, c, 4);CHKERRQ(ierr);
  }
  ierr = DMSetUp(*dm);CHKERRQ(ierr);
  for (c = 0; c < numCells; ++c) {
    GmfGetLin(meshIndex, GmfTetrahedra, &bufferTet[0], &bufferTet[1], &bufferTet[2], &bufferTet[3], &tag);
    for (i=0; i<4; ++i) bufferTet[i] += numCells - 1;  
    ierr = DMPlexSetCone(*dm, c, bufferTet);CHKERRQ(ierr);
  }

  ierr = DMSetDimension(*dm, dim);CHKERRQ(ierr);
  ierr = DMPlexSymmetrize(*dm);CHKERRQ(ierr);
  ierr = DMPlexStratify(*dm);CHKERRQ(ierr);

  
  DM idm = NULL;
  ierr = DMPlexInterpolate(*dm, &idm);CHKERRQ(ierr);
  ierr = DMDestroy(dm);CHKERRQ(ierr);
  *dm  = idm;

  if (bdLabelName) {
    ierr = PetscStrcmp(bdLabelName, "", &flg);CHKERRQ(ierr);
    if (!flg) {
      ierr = DMCreateLabel(*dm, bdLabelName);CHKERRQ(ierr);
      ierr = DMGetLabel(*dm, bdLabelName, &bdLabel);CHKERRQ(ierr);
    }
  }
  if (bdLabel) {
    numFacets = GmfStatKwd(meshIndex, GmfTriangles);
    GmfGotoKwd(meshIndex, GmfTriangles);
    for (f = 0; f < numFacets; ++f) {
      GmfGetLin(meshIndex, GmfTriangles, &bufferFac[0], &bufferFac[1], &bufferFac[2], &tag);
      for (i=0; i<3; ++i) bufferFac[i] += numCells - 1;
      ierr = DMPlexGetFullJoin(*dm, 3, (const PetscInt *) bufferFac, &joinSize, &join);CHKERRQ(ierr);
      if (joinSize != 1) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Could not determine Plex join for file triangle %d", f);
      ierr = DMLabelSetValue(bdLabel, join[0], tag);CHKERRQ(ierr);
      ierr = DMPlexRestoreJoin(*dm, 3, (const PetscInt *) bufferFac, &joinSize, &join);CHKERRQ(ierr);
    }

  }
  GmfCloseMesh(meshIndex);



  ierr = DMGetCoordinateSection(*dm, &coordSection);CHKERRQ(ierr);
  ierr = PetscSectionSetNumFields(coordSection, 1);CHKERRQ(ierr);
  ierr = PetscSectionSetFieldComponents(coordSection, 0, dim);CHKERRQ(ierr);
  ierr = PetscSectionSetChart(coordSection, numCells, numCells + numVertices);CHKERRQ(ierr);
  for (v = numCells; v < numCells+numVertices; ++v) {
    ierr = PetscSectionSetDof(coordSection, v, dim);CHKERRQ(ierr);
    ierr = PetscSectionSetFieldDof(coordSection, v, 0, dim);CHKERRQ(ierr);
  }
  ierr = PetscSectionSetUp(coordSection);CHKERRQ(ierr);
  ierr = PetscSectionGetStorageSize(coordSection, &coordSize);CHKERRQ(ierr);
  ierr = VecCreate(comm, &coordinates);CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject) coordinates, "coordinates");CHKERRQ(ierr);
  ierr = VecSetSizes(coordinates, coordSize, PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = VecSetType(coordinates, VECSTANDARD);CHKERRQ(ierr);
  ierr = VecGetArray(coordinates, &coords);CHKERRQ(ierr);
  for (v = 0; v < numVertices; ++v) {
    for (i = 0; i < dim; ++i) {
      coords[v*dim+i] = coordsIn[v*dim+i];
    }
  }
  ierr = VecRestoreArray(coordinates, &coords);CHKERRQ(ierr);
  ierr = PetscFree(coordsIn);CHKERRQ(ierr);
  ierr = DMSetCoordinatesLocal(*dm, coordinates);CHKERRQ(ierr);
  ierr = VecDestroy(&coordinates);CHKERRQ(ierr);


  PetscFunctionReturn(0);
}





#undef __FUNCT__
#define __FUNCT__ "DMPlexReadGmfSolFromFile_3d"
PetscErrorCode DMPlexReadGmfSolFromFile_3d(DM dm, PetscSection section, const char solName[], PetscInt solType, Vec * sol) {

  // TODO FOR NOW I ASSUME 1 SOL ONLY PER FILE
  
  MPI_Comm            comm;
  PetscScalar         * buffer;
  PetscInt            * ix;
  PetscInt            gmfVersion, dim, numSolAtVerticesLines, numSolTypes, solSize, locSolSize, solTypesTable[GmfMaxTyp];
  PetscInt            vStart, vEnd, numVertices, v, off, i;
  char                fileName[512];
  long long           solIndex;
  PetscErrorCode      ierr;


  PetscFunctionBegin;
  comm = PETSC_COMM_WORLD;
  ierr = DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd);CHKERRQ(ierr);
  numVertices = vEnd - vStart;
  ierr = VecCreate(comm, sol);CHKERRQ(ierr);
 
  strcpy(fileName, solName);
  strcat(fileName, ".solb");
  if ( !(solIndex = GmfOpenMesh(fileName, GmfRead, &gmfVersion, &dim)) ) {
    strcpy(fileName, solName);
    strcat(fileName,".sol");
    if ( !(solIndex = GmfOpenMesh(fileName, GmfRead, &gmfVersion, &dim)) ) {
      printf("####  ERROR Mesh file %s.mesh[b] not found ", solName);
      exit(1);
    }    
  }
  if (dim != 3) {
    printf("####  ERROR  Wrong dimension: %d != 3\n", dim);
    exit(1);
  }

  numSolAtVerticesLines = GmfStatKwd(solIndex, GmfSolAtVertices, &numSolTypes, &solSize, solTypesTable);  
  if( numSolAtVerticesLines == 0 ) {
    printf("####  ERROR  No SolAtVertices in the solution file %s\n", fileName);
    exit(1);
  }
  else if (numSolAtVerticesLines != numVertices) {
    printf("####  ERROR  The number of solution lines is different from the number of mesh vertices: %d != %d\n", numSolAtVerticesLines, numVertices);

  }
  if (numSolTypes > 1)
    printf("####  Warning  Several solution fields in file %s. Reading only the first one of type %d\n", fileName, solTypesTable[0]);
  
  if ((solTypesTable[0] != solType) && (solTypesTable[0] != 3 || solType != 5))
    printf("####  ERROR  Solution file solType and given solType not in agreement: %d != %d\n", solTypesTable[0], solType);  
  

  locSolSize = 0;
  switch (solType) {
    case 1: locSolSize = 1; break;
    case 2: locSolSize = 3; break;
    case 3: locSolSize = 9; break;
    case 5: locSolSize = 9; break;
    default: 
      printf("####  ERROR  non-scalar solutions not implemented yet\n");
      exit(1);  
  }
  ierr = VecSetSizes(*sol, PETSC_DECIDE, numSolAtVerticesLines*locSolSize);CHKERRQ(ierr);
  VecSetFromOptions(*sol);

  GmfGotoKwd(solIndex, GmfSolAtVertices);
  if (solType == 5)  {ierr = PetscMalloc2(solSize+3, &buffer, locSolSize, &ix);CHKERRQ(ierr);}
  else               {ierr = PetscMalloc2(solSize, &buffer, locSolSize, &ix);CHKERRQ(ierr);}
  for (v = 0; v < numSolAtVerticesLines; ++v) {
    GmfGetLin(solIndex, GmfSolAtVertices, buffer);
    ierr = PetscSectionGetOffset(section, v+vStart, &off);CHKERRQ(ierr);
    off /= dim;
    for (i=0; i<locSolSize; ++i) ix[i] = locSolSize*off + i;
    if (solType == 5) {
      // be careful, libmesh reads the matrix by columns
      buffer[8] = buffer[5]; buffer[7] = buffer[4]; buffer[6] = buffer[3]; 
      buffer[5] = buffer[4]; buffer[4] = buffer[2]; buffer[2] = buffer[3];
      buffer[3] = buffer[1];

    }
    VecSetValues(*sol, locSolSize, ix, buffer, INSERT_VALUES);
    
  }
  GmfCloseMesh(solIndex);
  ierr = PetscFree2(buffer, ix);CHKERRQ(ierr);
  VecAssemblyBegin(*sol);
  VecAssemblyEnd(*sol);

  PetscFunctionReturn(0);
}


