#include <petsc/private/dmfieldimpl.h> /*I "petscdmfield.h" I*/
#include <petsc/private/petscfeimpl.h> /*I "petscdmfield.h" I*/
#include <petscfe.h>
#include <petscdmplex.h>
#include <petscds.h>

typedef struct _n_DMField_DS
{
  PetscInt    fieldNum;
  Vec         vec;
  PetscInt    height;
  PetscObject *disc;
  PetscBool   multifieldVec;
}
DMField_DS;

static PetscErrorCode DMFieldDestroy_DS(DMField field)
{
  DMField_DS     *dsfield;
  PetscInt       i;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  dsfield = (DMField_DS *) field->data;
  ierr = VecDestroy(&dsfield->vec);CHKERRQ(ierr);
  for (i = 0; i < dsfield->height; i++) {
    ierr = PetscObjectDereference(dsfield->disc[i]);CHKERRQ(ierr);
  }
  ierr = PetscFree(dsfield->disc);CHKERRQ(ierr);
  ierr = PetscFree(dsfield);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode DMFieldView_DS(DMField field,PetscViewer viewer)
{
  DMField_DS     *dsfield = (DMField_DS *) field->data;
  PetscBool      iascii;
  PetscObject    disc;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscObjectTypeCompare((PetscObject)viewer,PETSCVIEWERASCII,&iascii);CHKERRQ(ierr);
  disc = dsfield->disc[0];
  if (iascii) {
    PetscViewerASCIIPrintf(viewer, "PetscDS field %D\n",dsfield->fieldNum);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
    ierr = PetscObjectView(disc,viewer);CHKERRQ(ierr);
    ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
  }
  ierr = PetscViewerASCIIPushTab(viewer);CHKERRQ(ierr);
  if (dsfield->multifieldVec) {
    SETERRQ(PetscObjectComm((PetscObject)field),PETSC_ERR_SUP,"View of subfield not implemented yet");
  } else {
    ierr = VecView(dsfield->vec,viewer);CHKERRQ(ierr);
  }
  ierr = PetscViewerASCIIPopTab(viewer);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode DMFieldEvaluate_DS(DMField field, Vec points, PetscDataType datatype, void *B, void *D, void *H)
{
  PetscFunctionBegin;
  SETERRQ(PetscObjectComm((PetscObject)field),PETSC_ERR_SUP,"Not implemented yet");
  PetscFunctionReturn(0);
}

static PetscErrorCode DMFieldDSGetHeightDisc(DMField field, PetscInt height, PetscObject *disc)
{
  DMField_DS     *dsfield = (DMField_DS *) field->data;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  if (!dsfield->disc[height]) {
    PetscClassId   id;

    ierr = PetscObjectGetClassId(dsfield->disc[0],&id);CHKERRQ(ierr);
    if (id == PETSCFE_CLASSID) {
      PetscFE fe = (PetscFE) dsfield->disc[0];

      ierr = PetscFECreateHeightTrace(fe,height,(PetscFE *)&dsfield->disc[height]);CHKERRQ(ierr);
    }
  }
  *disc = dsfield->disc[height];
  PetscFunctionReturn(0);
}

#define DMFieldDSdot(y,A,b,m,n,c,cast)                                           \
  do {                                                                           \
    PetscInt _i, _j, _k;                                                         \
    for (_i = 0; _i < (m); _i++) {                                               \
      for (_k = 0; _k < (c); _k++) {                                             \
        (y)[_i * (c) + _k] = 0.;                                                 \
      }                                                                          \
      for (_j = 0; _j < (n); _j++) {                                             \
        for (_k = 0; _k < (c); _k++) {                                           \
          (y)[_i * (c) + _k] += (A)[(_i * (n) + _j) * (c) + _k] * cast((b)[_j]); \
        }                                                                        \
      }                                                                          \
    }                                                                            \
  } while (0)

static PetscErrorCode DMFieldEvaluateFE_DS(DMField field, IS pointIS, PetscQuadrature quad, PetscDataType type, void *B, void *D, void *H)
{
  DMField_DS      *dsfield = (DMField_DS *) field->data;
  DM              dm;
  PetscObject     disc;
  PetscClassId    classid;
  PetscInt        nq, nc, dim, meshDim, numCells;
  PetscSection    section;
  const PetscReal *qpoints;
  PetscBool       isStride;
  const PetscInt  *points = NULL;
  PetscInt        sfirst = -1, stride = -1;
  PetscErrorCode  ierr;

  PetscFunctionBeginHot;
  dm   = field->dm;
  nc   = field->numComponents;
  ierr = PetscQuadratureGetData(quad,&dim,NULL,&nq,&qpoints,NULL);CHKERRQ(ierr);
  ierr = DMFieldDSGetHeightDisc(field,dsfield->height - 1 - dim,&disc);CHKERRQ(ierr);
  ierr = DMGetDimension(dm,&meshDim);CHKERRQ(ierr);
  ierr = DMGetDefaultSection(dm,&section);CHKERRQ(ierr);
  ierr = PetscSectionGetField(section,dsfield->fieldNum,&section);CHKERRQ(ierr);
  ierr = PetscObjectGetClassId(disc,&classid);CHKERRQ(ierr);
  /* TODO: batch */
  ierr = PetscObjectTypeCompare((PetscObject)pointIS,ISSTRIDE,&isStride);CHKERRQ(ierr);
  ierr = ISGetLocalSize(pointIS,&numCells);CHKERRQ(ierr);
  if (isStride) {
    ierr = ISStrideGetInfo(pointIS,&sfirst,&stride);CHKERRQ(ierr);
  } else {
    ierr = ISGetIndices(pointIS,&points);CHKERRQ(ierr);
  }
  if (classid == PETSCFE_CLASSID) {
    PetscFE      fe = (PetscFE) disc;
    PetscInt     feDim, i;
    PetscReal    *fB = NULL, *fD = NULL, *fH = NULL;

    if (dim == meshDim - 1) {
      /* TODO */
    }
    ierr = PetscFEGetDimension(fe,&feDim);CHKERRQ(ierr);
    ierr = PetscFEGetTabulation(fe,nq,qpoints,B ? &fB : NULL,D ? &fD : NULL,H ? &fH : NULL);CHKERRQ(ierr);
    for (i = 0; i < numCells; i++) {
      PetscInt     c = isStride ? (sfirst + i * stride) : points[i];
      PetscInt     closureSize;
      PetscScalar *elem = NULL;

      ierr = DMPlexVecGetClosure(dm,section,dsfield->vec,c,&closureSize,&elem);CHKERRQ(ierr);
      if (B) {
        if (type == PETSC_SCALAR) {
          PetscScalar *cB = &((PetscScalar *) B)[nc * nq * i];

          DMFieldDSdot(cB,fB,elem,nq,feDim,nc,(PetscScalar));
        } else {
          PetscReal *cB = &((PetscReal *) B)[nc * nq * i];

          DMFieldDSdot(cB,fB,elem,nq,feDim,nc,PetscRealPart);
        }
      }
      if (D) {
        if (type == PETSC_SCALAR) {
          PetscScalar *cD = &((PetscScalar *) D)[nc * nq * dim * i];

          DMFieldDSdot(cD,fD,elem,nq,feDim,(nc * dim),(PetscScalar));
        } else {
          PetscReal *cD = &((PetscReal *) D)[nc * nq * dim * i];

          DMFieldDSdot(cD,fD,elem,nq,feDim,(nc * dim),PetscRealPart);
        }
      }
      if (H) {
        if (type == PETSC_SCALAR) {
          PetscScalar *cH = &((PetscScalar *) H)[nc * nq * dim * dim * i];

          DMFieldDSdot(cH,fH,elem,nq,feDim,(nc * dim * dim),(PetscScalar));
        } else {
          PetscReal *cH = &((PetscReal *) H)[nc * nq * dim * dim * i];

          DMFieldDSdot(cH,fH,elem,nq,feDim,(nc * dim * dim),PetscRealPart);
        }
      }
      ierr = DMPlexVecRestoreClosure(dm,section,dsfield->vec,c,&closureSize,&elem);CHKERRQ(ierr);
    }
    ierr = PetscFERestoreTabulation(fe,nq,qpoints,B ? &fB : NULL,D ? &fD : NULL,H ? &fH : NULL);CHKERRQ(ierr);
  } else {SETERRQ(PetscObjectComm((PetscObject)field),PETSC_ERR_SUP,"Not implemented");}
  if (!isStride) {
    ierr = ISRestoreIndices(pointIS,&points);CHKERRQ(ierr);
  }
  ierr = PetscSectionDestroy(&section);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode DMFieldGetFEInvariance_DS(DMField field, IS pointIS, PetscBool *isConstant, PetscBool *isAffine, PetscBool *isQuadratic)
{
  DMField_DS     *dsfield;
  PetscObject    disc;
  PetscInt       h, imin;
  PetscClassId   id;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  dsfield = (DMField_DS *) field->data;
  ierr = ISGetMinMax(pointIS,&imin,NULL);CHKERRQ(ierr);
  for (h = 0; h < dsfield->height; h++) {
    PetscInt hEnd;

    ierr = DMPlexGetHeightStratum(field->dm,h,NULL,&hEnd);CHKERRQ(ierr);
    if (imin < hEnd) break;
  }
  ierr = DMFieldDSGetHeightDisc(field,h,&disc);CHKERRQ(ierr);
  ierr = PetscObjectGetClassId(disc,&id);CHKERRQ(ierr);
  if (id == PETSCFE_CLASSID) {
    PetscFE    fe = (PetscFE) disc;
    PetscInt   order, maxOrder;
    PetscBool  tensor = PETSC_FALSE;
    PetscSpace sp;

    ierr = PetscFEGetBasisSpace(fe, &sp);CHKERRQ(ierr);
    ierr = PetscSpaceGetOrder(sp,&order);CHKERRQ(ierr);
    ierr = PetscSpacePolynomialGetTensor(sp,&tensor);CHKERRQ(ierr);
    if (tensor) {
      PetscInt dim;

      ierr = DMGetDimension(field->dm,&dim);CHKERRQ(ierr);
      maxOrder = order * dim;
    } else {
      maxOrder = order;
    }
    if (isConstant)  *isConstant  = (maxOrder < 1) ? PETSC_TRUE : PETSC_FALSE;
    if (isAffine)    *isAffine    = (maxOrder < 2) ? PETSC_TRUE : PETSC_FALSE;
    if (isQuadratic) *isQuadratic = (maxOrder < 3) ? PETSC_TRUE : PETSC_FALSE;
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode DMFieldCreateDefaultQuadrature_DS(DMField field, IS pointIS, PetscQuadrature *quad)
{
  PetscInt       h, dim, imax, imin;
  DM             dm;
  DMField_DS     *dsfield;
  PetscObject    disc;
  PetscFE        fe;
  PetscClassId   id;
  PetscErrorCode ierr;


  PetscFunctionBegin;
  dm = field->dm;
  dsfield = (DMField_DS *) field->data;
  ierr = ISGetMinMax(pointIS,&imax,&imin);CHKERRQ(ierr);
  ierr = DMGetDimension(dm,&dim);CHKERRQ(ierr);
  for (h = 0; h <= dim; h++) {
    PetscInt hStart, hEnd;

    ierr = DMPlexGetHeightStratum(dm,h,&hStart,&hEnd);CHKERRQ(ierr);
    if (imin >= hStart && imax < hEnd) break;
  }
  *quad = NULL;
  if (h < dsfield->height) {
    ierr = DMFieldDSGetHeightDisc(field,h,&disc);CHKERRQ(ierr);
    ierr = PetscObjectGetClassId(disc,&id);CHKERRQ(ierr);
    if (id != PETSCFE_CLASSID) PetscFunctionReturn(0);
    fe = (PetscFE) disc;
    ierr = PetscFEGetQuadrature(fe,quad);CHKERRQ(ierr);
    ierr = PetscObjectReference((PetscObject)*quad);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

static PetscErrorCode DMFieldComputeFaceData_DS(DMField field, IS pointIS, PetscQuadrature quad, PetscFEGeom *geom)
{
  const PetscInt *points;
  PetscInt        p, dim, dE, numFaces, Nq;
  PetscBool       affineCells;
  DMLabel         depthLabel;
  IS              cellIS;
  DM              dm = field->dm;
  PetscErrorCode  ierr;

  PetscFunctionBegin;
  dim = geom->dim;
  dE  = geom->dimEmbed;
  ierr = DMPlexGetDepthLabel(dm, &depthLabel);CHKERRQ(ierr);
  ierr = DMLabelGetStratumIS(depthLabel, dim + 1, &cellIS);CHKERRQ(ierr);
  ierr = DMFieldGetFEInvariance(field,cellIS,NULL,&affineCells,NULL);CHKERRQ(ierr);
  ierr = ISGetIndices(pointIS, &points);CHKERRQ(ierr);
  numFaces = geom->numCells;
  Nq = geom->numPoints;
  if (affineCells) {
    PetscInt        numCells, offset, *cells;
    PetscFEGeom     *cellGeom;
    IS              suppIS;
    PetscQuadrature cellQuad = NULL;

    ierr = DMFieldCreateDefaultQuadrature(field,cellIS,&cellQuad);CHKERRQ(ierr);
    for (p = 0, numCells = 0; p < numFaces; p++) {
      PetscInt        point = points[p];
      PetscInt        numSupp, numChildren;

      ierr = DMPlexGetTreeChildren(dm, point, &numChildren, NULL); CHKERRQ(ierr);
      if (numChildren) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Face data not valid for facets with children");
      ierr = DMPlexGetSupportSize(dm, point,&numSupp);CHKERRQ(ierr);
      if (numSupp > 2) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Point %D has %D support, expected at most 2\n", point, numSupp);
      numCells += numSupp;
    }
    ierr = PetscMalloc1(numCells, &cells);CHKERRQ(ierr);
    for (p = 0, offset = 0; p < numFaces; p++) {
      PetscInt        point = points[p];
      PetscInt        numSupp, s;
      const PetscInt *supp;

      ierr = DMPlexGetSupportSize(dm, point,&numSupp);CHKERRQ(ierr);
      ierr = DMPlexGetSupport(dm, point, &supp);CHKERRQ(ierr);
      for (s = 0; s < numSupp; s++, offset++) {
        cells[offset] = supp[s];
      }
    }
    ierr = ISCreateGeneral(PETSC_COMM_SELF,numCells,cells,PETSC_USE_POINTER, &suppIS);CHKERRQ(ierr);
    ierr = DMFieldCreateFEGeom(field,suppIS,cellQuad,PETSC_FALSE,&cellGeom);CHKERRQ(ierr);
    for (p = 0, offset = 0; p < numFaces; p++) {
      PetscInt        point = points[p];
      PetscInt        numSupp, s, q;
      const PetscInt *supp;

      ierr = DMPlexGetSupportSize(dm, point,&numSupp);CHKERRQ(ierr);
      ierr = DMPlexGetSupport(dm, point, &supp);CHKERRQ(ierr);
      for (s = 0; s < numSupp; s++, offset++) {
        for (q = 0; q < Nq * dE * dE; q++) {
          geom->suppInvJ[s][p * Nq * dE * dE + q] = cellGeom->invJ[offset * Nq * dE * dE + q];
        }
      }
    }
    ierr = PetscFEGeomDestroy(&cellGeom);CHKERRQ(ierr);
    ierr = ISDestroy(&suppIS);CHKERRQ(ierr);
    ierr = PetscFree(cells);CHKERRQ(ierr);
    ierr = PetscQuadratureDestroy(&cellQuad);CHKERRQ(ierr);
  } else {
    PetscObject          faceDisc, cellDisc;
    PetscClassId         faceId, cellId;
    PetscDualSpace       dsp;
    DM                   K;
    PetscInt           (*co)[2][3];
    PetscInt             coneSize;
    PetscInt           **counts;
    PetscInt             f, i, o, q, s;
    const PetscInt      *coneK;
    PetscInt             minOrient, maxOrient, numOrient;
    PetscInt            *orients;
    PetscReal          **orientPoints;
    PetscReal           *cellPoints;
    PetscReal           *dummyWeights;
    PetscQuadrature      cellQuad = NULL;

    ierr = DMFieldDSGetHeightDisc(field, 1, &faceDisc);CHKERRQ(ierr);
    ierr = DMFieldDSGetHeightDisc(field, 0, &cellDisc);CHKERRQ(ierr);
    ierr = PetscObjectGetClassId(faceDisc,&faceId);CHKERRQ(ierr);
    ierr = PetscObjectGetClassId(cellDisc,&cellId);CHKERRQ(ierr);
    if (faceId != PETSCFE_CLASSID || cellId != PETSCFE_CLASSID) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Not supported\n");
    ierr = PetscFEGetDualSpace((PetscFE)cellDisc, &dsp);CHKERRQ(ierr);
    ierr = PetscDualSpaceGetDM(dsp, &K); CHKERRQ(ierr);
    ierr = DMPlexGetConeSize(K,0,&coneSize);CHKERRQ(ierr);
    ierr = DMPlexGetCone(K,0,&coneK);CHKERRQ(ierr);
    ierr = PetscMalloc4(numFaces,&co,dE*Nq,&cellPoints,coneSize,&counts,Nq,&dummyWeights);CHKERRQ(ierr);
    ierr = PetscQuadratureCreate(PetscObjectComm((PetscObject)field), &cellQuad);CHKERRQ(ierr);
    ierr = PetscQuadratureSetData(cellQuad, dE, 1, Nq, cellPoints, dummyWeights);CHKERRQ(ierr);
    minOrient = PETSC_MAX_INT;
    maxOrient = PETSC_MIN_INT;
    for (p = 0; p < numFaces; p++) { /* record the orientation of the facet wrt the support cells */
      PetscInt        point = points[p];
      PetscInt        numSupp, numChildren;
      const PetscInt *supp;

      ierr = DMPlexGetTreeChildren(dm, point, &numChildren, NULL); CHKERRQ(ierr);
      if (numChildren) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Face data not valid for facets with children");
      ierr = DMPlexGetSupportSize(dm, point,&numSupp);CHKERRQ(ierr);
      if (numSupp > 2) SETERRQ2(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Point %D has %D support, expected at most 2\n", point, numSupp);
      ierr = DMPlexGetSupport(dm, point, &supp);CHKERRQ(ierr);
      for (s = 0; s < numSupp; s++) {
        PetscInt        cell = supp[s];
        PetscInt        numCone;
        const PetscInt *cone, *orient;

        ierr = DMPlexGetConeSize(dm, cell, &numCone);CHKERRQ(ierr);
        if (numCone != coneSize) SETERRQ(PETSC_COMM_SELF, PETSC_ERR_PLIB, "Support point does not match reference element");
        ierr = DMPlexGetCone(dm, cell, &cone);CHKERRQ(ierr);
        ierr = DMPlexGetConeOrientation(dm, cell, &orient);CHKERRQ(ierr);
        for (f = 0; f < coneSize; f++) {
          if (cone[f] == point) break;
        }
        co[p][s][0] = f;
        co[p][s][1] = orient[f];
        co[p][s][2] = cell;
        minOrient = PetscMin(minOrient, orient[f]);
        maxOrient = PetscMin(maxOrient, orient[f]);
      }
      for (; s < 2; s++) {
        co[p][s][0] = -1;
        co[p][s][1] = -1;
        co[p][s][2] = -1;
      }
    }
    numOrient = maxOrient + 1 - minOrient;
    ierr = DMPlexGetCone(K,0,&coneK);CHKERRQ(ierr);
    /* count all (face,orientation) doubles that appear */
    ierr = PetscCalloc2(numOrient,&orients,numOrient,&orientPoints);CHKERRQ(ierr);
    for (f = 0; f < coneSize; f++) {ierr = PetscCalloc1(numOrient, &counts[f]);CHKERRQ(ierr);}
    for (p = 0; p < numFaces; p++) {
      for (s = 0; s < 2; s++) {
        if (co[p][s][0] >= 0) {
          counts[co[p][s][0]][co[p][s][1] - minOrient]++;
          orients[co[p][s][1] - minOrient]++;
        }
      }
    }
    for (o = 0; o < numOrient; o++) {
      if (orients[o]) {
        PetscInt orient = o + minOrient;
        PetscInt q;

        ierr = PetscMalloc1(Nq * dim, &orientPoints[o]);CHKERRQ(ierr);
        /* rotate the quadrature points appropriately */
        switch (dim) {
        case 0:
          break;
        case 1:
          if (orient == -2 || orient == 1) {
            for (q = 0; q < Nq; q++) {
              orientPoints[o][q] = -geom->xi[q];
            }
          } else {
            for (q = 0; q < Nq; q++) {
              orientPoints[o][q] = geom->xi[q];
            }
          }
          break;
        case 2:
          switch (coneSize) {
          case 3:
            for (q = 0; q < Nq; q++) {
              PetscReal lambda[3];
              PetscReal lambdao[3];

              /* convert to barycentric */
              lambda[0] = - (geom->xi[2 * q] + geom->xi[2 * q + 1]) / 2.;
              lambda[1] = (geom->xi[2 * q] + 1.) / 2.;
              lambda[2] = (geom->xi[2 * q + 1] + 1.) / 2.;
              if (orient >= 0) {
                for (i = 0; i < 3; i++) {
                  lambdao[i] = lambda[(orient + i) % 3];
                }
              } else {
                for (i = 0; i < 3; i++) {
                  lambdao[i] = lambda[(-(orient + i) + 3) % 3];
                }
              }
              /* convert to coordinates */
              orientPoints[o][2 * q + 0] = -(lambdao[0] + lambdao[2]) + lambdao[1];
              orientPoints[o][2 * q + 1] = -(lambdao[0] + lambdao[1]) + lambdao[2];
            }
            break;
          case 4:
            for (q = 0; q < Nq; q++) {
              PetscReal xi[2], xio[2];
              PetscInt oabs = (orient >= 0) ? orient : -(orient + 1);

              xi[0] = geom->xi[2 * q];
              xi[1] = geom->xi[2 * q + 1];
              switch (oabs) {
              case 0:
                xio[0] = xi[0];
                xio[1] = xi[1];
                break;
              case 1:
                xio[0] = xi[1];
                xio[1] = -xi[0];
                break;
              case 2:
                xio[0] = -xi[0];
                xio[1] = -xi[1];
              case 3:
                xio[0] = -xi[1];
                xio[1] = xi[0];
              }
              if (orient < 0) {
                xio[0] = -xio[0];
              }
              orientPoints[o][2 * q + 0] = xio[0];
              orientPoints[o][2 * q + 1] = xio[1];
            }
            break;
          default:
            SETERRQ(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Not implemented yet\n");
          }
        default:
          SETERRQ(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Not implemented yet\n");
        }
      }
    }
    for (f = 0; f < coneSize; f++) {
      PetscInt face = coneK[f];
      PetscScalar v0[3];
      PetscScalar J[9];
      PetscInt numCells, offset;
      PetscInt *cells;
      IS suppIS;

      ierr = DMPlexComputeCellGeometryFEM(K, face, NULL, v0, J, NULL, NULL);CHKERRQ(ierr);
      for (o = 0; o <= numOrient; o++) {
        PetscFEGeom *cellGeom;

        if (!counts[f][o]) continue;
        /* If this (face,orientation) double appears,
         * convert the face quadrature points into volume quadrature points */
        for (q = 0; q < Nq; q++) {
          PetscReal xi0[3] = {-1., -1., -1.};

          CoordinatesRefToReal(dE, dim, xi0, v0, J, &orientPoints[o][dim * q + 0], &cellPoints[dE * q + 0]);
        }
        for (p = 0, numCells = 0; p < numFaces; p++) {
          for (s = 0; s < 2; s++) {
            if (co[p][s][0] == f && co[p][s][1] == o + minOrient) numCells++;
          }
        }
        ierr = PetscMalloc1(numCells, &cells);CHKERRQ(ierr);
        for (p = 0, offset = 0; p < numFaces; p++) {
          for (s = 0; s < 2; s++) {
            if (co[p][s][0] == f && co[p][s][1] == o + minOrient) {
              cells[offset++] = co[p][s][2];
            }
          }
        }
        ierr = ISCreateGeneral(PETSC_COMM_SELF,numCells,cells,PETSC_USE_POINTER, &suppIS);CHKERRQ(ierr);
        ierr = DMFieldCreateFEGeom(field,suppIS,cellQuad,PETSC_FALSE,&cellGeom);CHKERRQ(ierr);
        for (p = 0, offset = 0; p < numFaces; p++) {
          for (s = 0; s < 2; s++) {
            if (co[p][s][0] == f && co[p][s][1] == o + minOrient) {
              for (q = 0; q < Nq * dE * dE; q++) {
                geom->suppInvJ[s][p * Nq * dE * dE + q] = cellGeom->invJ[offset * Nq * dE * dE + q];
              }
              offset++;
            }
          }
        }
        ierr = PetscFEGeomDestroy(&cellGeom);CHKERRQ(ierr);
        ierr = ISDestroy(&suppIS);CHKERRQ(ierr);
        ierr = PetscFree(cells);CHKERRQ(ierr);
      }
    }
    for (o = 0; o < numOrient; o++) {
      if (orients[o]) {
        ierr = PetscFree(orientPoints[o]);CHKERRQ(ierr);
      }
    }
    ierr = PetscFree2(orients,orientPoints);CHKERRQ(ierr);
    ierr = PetscQuadratureDestroy(&cellQuad);CHKERRQ(ierr);
    ierr = PetscFree4(co,cellPoints,counts,dummyWeights);CHKERRQ(ierr);
  }
  ierr = ISRestoreIndices(pointIS, &points);CHKERRQ(ierr);
  ierr = ISDestroy(&cellIS);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

static PetscErrorCode DMFieldInitialize_DS(DMField field)
{
  PetscFunctionBegin;
  field->ops->destroy                 = DMFieldDestroy_DS;
  field->ops->evaluate                = DMFieldEvaluate_DS;
  field->ops->evaluateFE              = DMFieldEvaluateFE_DS;
  field->ops->getFEInvariance         = DMFieldGetFEInvariance_DS;
  field->ops->createDefaultQuadrature = DMFieldCreateDefaultQuadrature_DS;
  field->ops->view                    = DMFieldView_DS;
  field->ops->computeFaceData         = DMFieldComputeFaceData_DS;
  PetscFunctionReturn(0);
}

PETSC_INTERN PetscErrorCode DMFieldCreate_DS(DMField field)
{
  DMField_DS     *dsfield;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscNewLog(field,&dsfield);CHKERRQ(ierr);
  field->data = dsfield;
  ierr = DMFieldInitialize_DS(field);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

PetscErrorCode DMFieldCreateDS(DM dm, PetscInt fieldNum, Vec vec,DMField *field)
{
  DMField        b;
  DMField_DS     *dsfield;
  PetscObject    disc = NULL;
  PetscBool      isContainer = PETSC_FALSE;
  PetscClassId   id = -1;
  PetscInt       numComponents = -1, dsNumFields;
  PetscSection   section;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = DMGetDefaultSection(dm,&section);CHKERRQ(ierr);
  ierr = PetscSectionGetFieldComponents(section,fieldNum,&numComponents);CHKERRQ(ierr);
  ierr = DMGetNumFields(dm,&dsNumFields);CHKERRQ(ierr);
  if (dsNumFields) {ierr = DMGetField(dm,fieldNum,&disc);CHKERRQ(ierr);}
  if (disc) {
    ierr = PetscObjectGetClassId(disc,&id);CHKERRQ(ierr);
    isContainer = (id == PETSC_CONTAINER_CLASSID) ? PETSC_TRUE : PETSC_FALSE;
  }
  if (!disc || isContainer) {
    MPI_Comm        comm = PetscObjectComm((PetscObject) dm);
    PetscInt        cStart, cEnd, dim;
    PetscInt        localConeSize = 0, coneSize;
    PetscFE         fe;
    PetscDualSpace  Q;
    PetscSpace      P;
    DM              K;
    PetscQuadrature quad, fquad;
    PetscBool       isSimplex;

    ierr = DMPlexGetHeightStratum(dm, 0, &cStart, &cEnd);CHKERRQ(ierr);
    ierr = DMGetDimension(dm, &dim);CHKERRQ(ierr);
    if (cEnd > cStart) {
      ierr = DMPlexGetConeSize(dm, cStart, &localConeSize);CHKERRQ(ierr);
    }
    ierr = MPI_Allreduce(&localConeSize,&coneSize,1,MPIU_INT,MPI_MAX,comm);CHKERRQ(ierr);
    isSimplex = (coneSize == (dim + 1)) ? PETSC_TRUE : PETSC_FALSE;
    ierr = PetscSpaceCreate(comm, &P);CHKERRQ(ierr);
    ierr = PetscSpaceSetOrder(P, 1);CHKERRQ(ierr);
    ierr = PetscSpaceSetNumComponents(P, numComponents);CHKERRQ(ierr);
    ierr = PetscSpaceSetType(P,PETSCSPACEPOLYNOMIAL);CHKERRQ(ierr);
    ierr = PetscSpaceSetNumVariables(P, dim);CHKERRQ(ierr);
    ierr = PetscSpacePolynomialSetTensor(P, isSimplex ? PETSC_FALSE : PETSC_TRUE);CHKERRQ(ierr);
    ierr = PetscSpaceSetUp(P);CHKERRQ(ierr);
    ierr = PetscDualSpaceCreate(comm, &Q);CHKERRQ(ierr);
    ierr = PetscDualSpaceSetType(Q,PETSCDUALSPACELAGRANGE);CHKERRQ(ierr);
    ierr = PetscDualSpaceCreateReferenceCell(Q, dim, isSimplex, &K);CHKERRQ(ierr);
    ierr = PetscDualSpaceSetDM(Q, K);CHKERRQ(ierr);
    ierr = DMDestroy(&K);CHKERRQ(ierr);
    ierr = PetscDualSpaceSetNumComponents(Q, numComponents);CHKERRQ(ierr);
    ierr = PetscDualSpaceSetOrder(Q, 1);CHKERRQ(ierr);
    ierr = PetscDualSpaceLagrangeSetTensor(Q, isSimplex ? PETSC_FALSE : PETSC_TRUE);CHKERRQ(ierr);
    ierr = PetscDualSpaceSetUp(Q);CHKERRQ(ierr);
    ierr = PetscFECreate(comm, &fe);CHKERRQ(ierr);
    ierr = PetscFESetType(fe,PETSCFEBASIC);CHKERRQ(ierr);
    ierr = PetscFESetBasisSpace(fe, P);CHKERRQ(ierr);
    ierr = PetscFESetDualSpace(fe, Q);CHKERRQ(ierr);
    ierr = PetscFESetNumComponents(fe, numComponents);CHKERRQ(ierr);
    ierr = PetscFESetUp(fe);CHKERRQ(ierr);
    ierr = PetscSpaceDestroy(&P);CHKERRQ(ierr);
    ierr = PetscDualSpaceDestroy(&Q);CHKERRQ(ierr);
    if (isSimplex) {
      ierr = PetscDTGaussJacobiQuadrature(dim,   1, 1, -1.0, 1.0, &quad);CHKERRQ(ierr);
      ierr = PetscDTGaussJacobiQuadrature(dim-1, 1, 1, -1.0, 1.0, &fquad);CHKERRQ(ierr);
    }
    else {
      ierr = PetscDTGaussTensorQuadrature(dim,   1, 1, -1.0, 1.0, &quad);CHKERRQ(ierr);
      ierr = PetscDTGaussTensorQuadrature(dim-1, 1, 1, -1.0, 1.0, &fquad);CHKERRQ(ierr);
    }
    ierr = PetscFESetQuadrature(fe, quad);CHKERRQ(ierr);
    ierr = PetscFESetFaceQuadrature(fe, fquad);CHKERRQ(ierr);
    ierr = PetscQuadratureDestroy(&quad);CHKERRQ(ierr);
    ierr = PetscQuadratureDestroy(&fquad);CHKERRQ(ierr);
    disc = (PetscObject) fe;
  } else {
    ierr = PetscObjectReference(disc);CHKERRQ(ierr);
  }
  ierr = PetscObjectGetClassId(disc,&id);CHKERRQ(ierr);
  if (id == PETSCFE_CLASSID) {
    PetscFE fe = (PetscFE) disc;

    ierr = PetscFEGetNumComponents(fe,&numComponents);CHKERRQ(ierr);
  } else {SETERRQ(PetscObjectComm((PetscObject)dm),PETSC_ERR_SUP,"Not implemented");}
  ierr = DMFieldCreate(dm,numComponents,DMFIELD_VERTEX,&b);CHKERRQ(ierr);
  ierr = DMFieldSetType(b,DMFIELDDS);CHKERRQ(ierr);
  dsfield = (DMField_DS *) b->data;
  dsfield->fieldNum = fieldNum;
  ierr = DMGetDimension(dm,&dsfield->height);CHKERRQ(ierr);
  dsfield->height++;
  ierr = PetscCalloc1(dsfield->height,&dsfield->disc);CHKERRQ(ierr);
  dsfield->disc[0] = disc;
  ierr = PetscObjectReference((PetscObject)vec);CHKERRQ(ierr);
  dsfield->vec = vec;
  *field = b;

  PetscFunctionReturn(0);
}
