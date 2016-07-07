#include <petsc/private/dmpleximpl.h>   /*I      "petscdmplex.h"   I*/
#include "libmeshb7.h"




#undef __FUNCT__
#define __FUNCT__ "DMPlexWrite_gmfMesh2d_1sol"
// Temporary - For simple python use
PetscErrorCode DMPlexWrite_gmfMesh2d_1sol(DM dm, PetscBool writeMesh, const char bdLabelName[], const char meshName[], 
                                          Vec sol,  PetscInt solType, const char solNames[], 
                                          PetscSection section, PetscBool ascii) {
  
  PetscErrorCode ierr;
        
  PetscFunctionBegin;
  ierr = DMPlexWrite_gmfMesh2d(dm, writeMesh, bdLabelName, meshName, 1, &sol,  &solType, &solNames, section, ascii);CHKERRQ(ierr);                                
  PetscFunctionReturn(0);                                      

}


#undef __FUNCT__
#define __FUNCT__ "DMPlexWrite_gmfMesh2d_noSol"
// Temporary - For simple python use
PetscErrorCode DMPlexWrite_gmfMesh2d_noSol(DM dm, const char bdLabelName[], const char meshName[], PetscSection section, PetscBool ascii) {
  
  PetscErrorCode ierr;
        
  PetscFunctionBegin;
  ierr = DMPlexWrite_gmfMesh2d(dm, PETSC_TRUE, bdLabelName, meshName, 0, NULL, NULL, NULL, section, ascii);CHKERRQ(ierr);                                
  PetscFunctionReturn(0);                                      

}


#undef __FUNCT__
#define __FUNCT__ "DMPlexWrite_gmfMesh2d"
// solTypes: 1 for a scalar, 2 for a vector, 3 for symmetric matrix and 4 for a full matrix

PetscErrorCode DMPlexWrite_gmfMesh2d(DM dm, PetscBool writeMesh, const char bdLabelName[], const char meshName[], 
                                     PetscInt numSol, Vec * sol,  PetscInt * solTypes, const char * solNames[],
                                     PetscSection section, PetscBool ascii) {
  
  DM                  cdm;
  PetscSection        coordSection;
  Vec                 coordinates;
  DMLabel             bdLabel = NULL;
  const PetscScalar * coords, * solution;
  PetscScalar       * buffer;
  PetscInt            dim, cStart, cEnd, numCells, c, vStart, vEnd, numVertices, v, eStart, eEnd, numEdges, e, off;
  PetscInt            idx[3], i, iSol, tag, coneSize;
  const PetscInt    * cone;
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
  ierr = DMPlexGetDepthStratum(dm, 1, &eStart, &eEnd);CHKERRQ(ierr);
  numEdges = eEnd - eStart;
  if (section)
    coordSection = section;
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
      GmfSetLin(meshIndex, GmfVertices, coords[off], coords[off+1], tag);  
    }
    ierr = VecRestoreArrayRead(coordinates, &coords);CHKERRQ(ierr);
    
    GmfSetKwd(meshIndex, GmfTriangles, numCells);
    tag = 0;
    for (c = cStart; c < cEnd; ++c) {
      PetscInt *closure = NULL;
      PetscInt  closureSize, cl, p;
      ierr = DMPlexGetTransitiveClosure(dm, c, PETSC_TRUE, &closureSize, &closure);CHKERRQ(ierr);
      for (cl = 0, i=0; cl < closureSize*2; cl += 2) {
        p = closure[cl];
        if (p >= vStart && p < vEnd) 
          idx[i++] = p-vStart+1;
      }
      GmfSetLin(meshIndex, GmfTriangles, idx[0], idx[1], idx[2], tag);  
      ierr = DMPlexRestoreTransitiveClosure(dm, c, PETSC_TRUE, &closureSize, &closure);CHKERRQ(ierr);
    }

    GmfSetKwd(meshIndex, GmfEdges, numEdges);
    if (bdLabelName) {
      ierr = PetscStrcmp(bdLabelName, "", &flg);CHKERRQ(ierr);
      if (!flg) {
        ierr = DMGetLabel(dm, bdLabelName, &bdLabel);CHKERRQ(ierr);
      }
    }
    if (bdLabel)  {
      for (e = eStart; e < eEnd; ++e) {
        ierr = DMPlexGetConeSize(dm, e, &coneSize);CHKERRQ(ierr);
        ierr = DMPlexGetCone(dm, e, &cone);CHKERRQ(ierr);
        if (coneSize != 2) {
          printf("####  ERROR  cone of an Edge != 2, edge: %d, size: %d\n", e, coneSize);
          exit(1);
        }
        idx[0] = cone[0] - vStart + 1; idx[1] = cone[1] - vStart + 1;
        ierr = DMLabelGetValue(bdLabel, e, &tag);CHKERRQ(ierr);
        if (tag>0) GmfSetLin(meshIndex, GmfEdges, idx[0], idx[1], tag);
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

    GmfSetKwd(solIndex, solKeyword, numVertices, 1, &solTypes[iSol]);
    ierr = VecGetArrayRead(sol[iSol], &solution);CHKERRQ(ierr);
    buffer = NULL;
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
    if ( !GmfCloseMesh(solIndex) ) {
      fprintf(stderr,"#### ERROR: mesh file %s cannot be closed\n",fileName);
      exit(1);
    }

  }
  


  PetscFunctionReturn(0);
  
}








#undef __FUNCT__
#define __FUNCT__ "DMPlexCreateGmfFromFile_2d"
PetscErrorCode DMPlexCreateGmfFromFile_2d(const char meshName[], const char bdLabelName[], DM *dm)
{

  MPI_Comm            comm;
  PetscSection        coordSection;
  Vec                 coordinates;
  DMLabel             bdLabel = NULL;
  PetscScalar         * coordsIn, * coords;
  PetscScalar         buffer[2];
  PetscInt            bufferTri[3], bufferEdg[2];
  PetscInt            dim, numCells, c, numVertices, v, numEdges, e, coordSize;
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
  if (dim != 2) {
    printf("####  ERROR  Wrong dimension: %d != 2\n", dim);
    exit(1);
  }

  numVertices = GmfStatKwd(meshIndex, GmfVertices);
  numCells    = GmfStatKwd(meshIndex, GmfTriangles);

  if (numVertices <= 0 ) {
    printf("####  ERROR  Number of vertices: %d <= 0\n", numVertices);
    exit(1);
  }

  ierr = PetscMalloc1(numVertices*dim, &coordsIn);CHKERRQ(ierr);
  GmfGotoKwd(meshIndex, GmfVertices);
  for (v = 0; v < numVertices ; ++v) {
    //printf("DEBUG  v: %d\n", v);
    GmfGetLin(meshIndex, GmfVertices, &buffer[0], &buffer[1], &tag);
    coordsIn[2*v]    = (double)buffer[0];
    coordsIn[2*v+1]  = (double)buffer[1];
  }

  GmfGotoKwd(meshIndex, GmfTriangles);
  ierr = DMPlexSetChart(*dm, 0, numCells+numVertices);CHKERRQ(ierr);
  for (c = 0; c < numCells; ++c) {
    ierr = DMPlexSetConeSize(*dm, c, 3);CHKERRQ(ierr);
  }
  ierr = DMSetUp(*dm);CHKERRQ(ierr);
  for (c = 0; c < numCells; ++c) {
    GmfGetLin(meshIndex, GmfTriangles, &bufferTri[0], &bufferTri[1], &bufferTri[2], &tag);
    for (i=0; i<3; ++i) bufferTri[i] += numCells - 1;  
    ierr = DMPlexSetCone(*dm, c, bufferTri);CHKERRQ(ierr);
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
    numEdges    = GmfStatKwd(meshIndex, GmfEdges);
    GmfGotoKwd(meshIndex, GmfEdges);
    for (e = 0; e < numEdges; ++e) {
      GmfGetLin(meshIndex, GmfEdges, &bufferEdg[0], &bufferEdg[1], &tag);
      bufferEdg[0] += numCells - 1; bufferEdg[1] += numCells - 1;
      ierr = DMPlexGetFullJoin(*dm, 2, (const PetscInt *) bufferEdg, &joinSize, &join);CHKERRQ(ierr);
      if (joinSize != 1) SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "Could not determine Plex join for file edge %d", e);
      ierr = DMLabelSetValue(bdLabel, join[0], tag);CHKERRQ(ierr);
      ierr = DMPlexRestoreJoin(*dm, 2, (const PetscInt *) bufferEdg, &joinSize, &join);CHKERRQ(ierr);
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
#define __FUNCT__ "DMPlexReadGmfSolFromFile_2d"
PetscErrorCode DMPlexReadGmfSolFromFile_2d(DM dm, PetscSection section, const char solName[], Vec * sol, PetscInt * solType) {

  // TODO FOR NOW I ASSUME 1 SOL ONLY PER FILE
  
  MPI_Comm            comm;
  PetscScalar         * buffer;
  PetscInt            * ix;
  PetscInt            gmfVersion, dim, numSolAtVerticesLines, numSolTypes, solSize, solTypesTable[GmfMaxTyp];
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
  if (dim != 2) {
    printf("####  ERROR  Wrong dimension: %d != 2\n", dim);
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
  
  
  VecSetSizes(*sol, PETSC_DECIDE, numSolAtVerticesLines*solTypesTable[0]);
  VecSetFromOptions(*sol);


  GmfGotoKwd(solIndex, GmfSolAtVertices);
  ierr = PetscMalloc2(solSize, &buffer, solTypesTable[0], &ix);CHKERRQ(ierr);
  for (v = 0; v < numSolAtVerticesLines; ++v) {
    GmfGetLin(solIndex, GmfSolAtVertices, buffer);
    ierr = PetscSectionGetOffset(section, v+vStart, &off);CHKERRQ(ierr);
    for (i=0; i<solTypesTable[0]; ++i) ix[i] = solTypesTable[0]*off/2 + i;
    VecSetValues(*sol, solTypesTable[0], ix, buffer, INSERT_VALUES);
  }
  ierr = PetscFree2(buffer, ix);CHKERRQ(ierr);
  VecAssemblyBegin(*sol);
  VecAssemblyEnd(*sol);

  *solType = solTypesTable[0];

  PetscFunctionReturn(0);
}


