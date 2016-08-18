#include <petsc/private/dmpleximpl.h>   /*I      "petscdmplex.h"   I*/



#undef __FUNCT__
#define __FUNCT__ "DMPlexMetricReduction2d_Internal"
// metric given as upper triangular matrix mat
PetscErrorCode DMPlexMetricReduction2d_Internal(PetscReal * mat, PetscReal * eigVal, PetscReal * eigVec) {
  
  PetscInt    i;
  PetscReal   nrm, sign, vecNrm, dd, tmp;
  PetscReal met[3];
  
  PetscFunctionBegin;
  
  met[0] = mat[0]; met[1] = mat[1]; met[2] = mat[2];
    
  nrm = fabs(met[0]);
  if ( met[0] >= 0. ) sign =  1.;
  else sign = -1;  
  for (i=1; i<3; ++i) {
    if ( fabs(met[i]) > nrm ) {
      nrm = fabs(met[i]);
      if ( met[i] >= 0. ) sign =  1.;
      else  sign = -1;  
    }  
  }    
    
  //--- a null matrix
  if ( nrm < 1e-100 ) {
    eigVal[0] = eigVal[1]= 0.;
    eigVec[0] = 1.; eigVec[1] = 0.;
    eigVec[2] = 0.; eigVec[3] = 1.;
    PetscFunctionReturn(0);
  }
  
  nrm = sign*nrm;
  
  //--- compute eigenvalues
  for (i=0; i<3; ++i) met[i] = met[i]/nrm;
  
  dd = (met[0]-met[2])*(met[0]-met[2]) + 4.*met[1]*met[1];

  //--- Diagonal matrix with two identical eigenvalues
  if ( fabs(dd) < 1.e-24 ) {  // eigVal[0] = eigVal[1] = nrm (or 1 after normalization)
    eigVal[0] = eigVal[1] = nrm;
    eigVec[0] = 1.; eigVec[1] = 0.;
    eigVec[2] = 0.; eigVec[3] = 1.;
    PetscFunctionReturn(0);
  }
  else if ( dd < 0. ) exit(10);
  
  dd = sqrt(dd);
  eigVal[0] = 0.5*(met[0]+met[2]-dd);
  eigVal[1] = 0.5*(met[0]+met[2]+dd);
  
  //--- compute eigenvectors
  eigVec[2] = -met[1];
  eigVec[3] =  met[0]-eigVal[1];
  vecNrm    =  eigVec[2]*eigVec[2] + eigVec[3]*eigVec[3];

  if ( vecNrm < 1e-30 ) { // => diag matrix => bad sol use the other
    eigVec[2] = eigVal[1]-met[2];
    eigVec[3] = met[1];
    vecNrm    = eigVec[2]*eigVec[2] + eigVec[3]*eigVec[3];
    if ( vecNrm < 1e-30 ) { //--- in the case we have dd = 0 ie egv1 = egv2 
    //--- thus M is the Id matrix after normalisation
      if ( fabs(eigVal[0]-1.) < 1.e-12 && fabs(eigVal[1]-1.) < 1.e-12 ) {  
      eigVal[0] = eigVal[1] = nrm;
        eigVec[0] = 1.; eigVec[1] = 0.;
        eigVec[2] = 0.; eigVec[3] = 1.;
      PetscFunctionReturn(0);
      }
      else
        exit(13);
    }
  }
  
  vecNrm     = 1./sqrt(vecNrm);
  eigVec[2] *= vecNrm;
  eigVec[3] *= vecNrm;

  if ( fabs(eigVal[0]) < fabs(eigVal[1]) ) {
    eigVec[0] =  eigVec[3];
    eigVec[1] = -eigVec[2];
  }
  else {
    eigVec[0] =  eigVec[2];
    eigVec[1] =  eigVec[3];
    eigVec[2] = -eigVec[1];
    eigVec[3] =  eigVec[0];
    tmp       = eigVal[0];
    eigVal[0] = eigVal[1];
    eigVal[1] = tmp;
  }

  eigVal[0] *= nrm;
  eigVal[1] *= nrm;

  PetscFunctionReturn(0);
  
}

static inline void intersectDegeneratedCase2D(double *eigVal1, double *eigVec1, double *eigVal2, double *eigVec2,
                                                  double *metNew)
{
  int    pos1, pos2;       // position of the non null e.v. in each metric
  int    posNul1, posNul2; // position of the null e.v. in each metric
  double val,ps1,ps2;
  double eigVal[2],eigVec[4],detPinv,Pinv[4];
  
  double eps = 1e-40;


  pos1 = (fabs(eigVal1[0]) > eps ) ? 0 : 1;
  pos2 = (fabs(eigVal2[0]) > eps ) ? 0 : 1;
  
  //--- if det = 0, ie if the two eigen vectors corresponding to the non null eigenvalues are colinear
  val = eigVec1[2*pos1+0]*eigVec2[2*pos2+1] - eigVec1[2*pos1+1]*eigVec2[2*pos2+0];
  
  if ( val < eps*eps ) {
    eigVec[0] =  eigVec1[2*pos1+0];
    eigVec[1] =  eigVec1[2*pos1+1];
    eigVec[2] = -eigVec1[2*pos1+1];  
    eigVec[3] =  eigVec1[2*pos1+0];
    eigVal[0] =  fmax(eigVal1[pos1], eigVal2[pos2]);
    eigVal[1] =  0;
  }
  else {
    // position of the null eigenvalues:
    posNul1   = abs(pos1-1);
    posNul2   = abs(pos2-1);
    eigVec[0] = eigVec1[2*posNul1+0];
    eigVec[1] = eigVec1[2*posNul1+1];
    eigVec[2] = eigVec2[2*posNul2+0];   
    eigVec[3] = eigVec2[2*posNul2+1];
    ps1       = eigVec1[2*posNul1+0]*eigVec2[2*pos2+0] + eigVec1[2*posNul1+1]*eigVec2[2*pos2+1];
    ps2       = eigVec1[2*pos1+0]*eigVec2[2*posNul2+0] + eigVec1[2*pos1+1]*eigVec2[2*posNul2+1];
    eigVal[1] = eigVal1[pos1]*ps1*ps1;
    eigVal[0] = eigVal2[pos2]*ps2*ps2;
  }

  //--- Inversion of the transformation matrix
  detPinv =  1./(eigVec[0]*eigVec[3]-eigVec[1]*eigVec[2]);
  Pinv[0] =  eigVec[3]*detPinv;
  Pinv[1] = -eigVec[2]*detPinv;
  Pinv[2] = -eigVec[1]*detPinv;
  Pinv[3] =  eigVec[0]*detPinv;

  metNew[0] = eigVal[0]*Pinv[0]*Pinv[0] + eigVal[1]*Pinv[2]*Pinv[2];
  metNew[1] = eigVal[0]*Pinv[0]*Pinv[1] + eigVal[1]*Pinv[2]*Pinv[3];
  metNew[2] = eigVal[0]*Pinv[1]*Pinv[1] + eigVal[1]*Pinv[3]*Pinv[3];
  
  return;
}



PetscErrorCode metricIntersection(double* metric1, double* metric2, double* metNew)
{
  int    i,cpt[2];
  double met1[3],met2[3];
  double eigVal[2], eigVec[4] = {0,0,0,0};
  double tmp;
  double eigVal1[2] = {0,0}, eigVal2[2] = {0,0};
  double eigVec1[4] = {0,0,0,0}, eigVec2[4] = {0,0,0,0};
  double sqrtMat[3], sqrtMatInv[3], M2bar[3], Mintbar[3]; 
  
  PetscInt ierr;


  double eps = 1e-40;

  PetscFunctionBegin;

  //--- Copy metric in working array because they may be switched
  for (i=0; i<3; ++i) {    
    met1[i] = metric1[i]; 
    met2[i] = metric2[i]; 
  }

  //--- compute the eigen vectors of N=M1^(-1)*M2
  //--- we have: vp such that : M1v = vp*M2v <=> (M1-vp*M2)v = 0 
  ierr = DMPlexMetricReduction2d_Internal(met1, eigVal1, eigVec1);CHKERRQ(ierr);
  ierr = DMPlexMetricReduction2d_Internal(met2, eigVal2, eigVec2);CHKERRQ(ierr);

  // check for degenerated cases
  // -0- Normal cases (M1 and M2 are > 0 ) 
  // -1- M1 and M2 are zeros                => metNew = 0
  // -2- If M1 is zero                      => metNew = M2
  // -3- If M2 is zero                      => metNew = M1
  // -4- If M1 has 1 zero eigVal and not M2 => switch M1 and M2
  // -5- If M2 has 1 zero eigVal and not M1 => nothing to do
  // -6- If each matrix has a zero eigVal   => specific treatment
  
  cpt[0] = cpt[1] = 0;
  for (i=0; i<2; ++i) {
    if ( fabs(eigVal1[i]) < eps ) cpt[0]++;
    if ( fabs(eigVal2[i]) < eps ) cpt[1]++;
  }
  
  if ( cpt[0] == 0 ) {
    if ( cpt[1] != 2 ) {
      // do nothing, normal case
      for (i=0; i<2; ++i) eigVal[i] = eigVal1[i];
      for (i=0; i<4; ++i) eigVec[i] = eigVec1[i];
    }
    else {
      // met2 is null
      for (i=0; i<3; ++i) metNew[i] = met1[i];
      PetscFunctionReturn(0);           
    }
  }
  else if ( cpt[0] == 1 ) {
    if ( cpt[1] == 0 ) {
      // switch
      for (i=0; i<3; ++i){
        tmp     = met1[i];
        met1[i] = met2[i];
        met2[i] = tmp;
      }
      for (i=0; i<2; ++i)
        eigVal[i] = eigVal2[i];
      for (i=0; i<4; ++i)
        eigVec[i] = eigVec2[i];
    }
    else if ( cpt[1] == 1 ) {
      //--- Treat specific degenerated case
      intersectDegeneratedCase2D(eigVal1,eigVec1,eigVal2,eigVec2,metNew);
      PetscFunctionReturn(0);
    }
    else {
      for (i=0; i<3; ++i) metNew[i] = met1[i];
      PetscFunctionReturn(0);           
    }
  }
  else {
    if ( cpt[1] != 2 ) {
      for (i=0; i<3; ++i) metNew[i] = met2[i];
      PetscFunctionReturn(0);           
    }
    else {
      for (i=0; i<3; ++i)  metNew[i] = 0.;
      PetscFunctionReturn(0);           
    }
  }
  
  
  //--- Compute sqrt(Lbd_M1)
  eigVal[0] = sqrt(fabs(eigVal[0]));
  eigVal[1] = sqrt(fabs(eigVal[1]));  

  //--- Compute M1^(1/2)  (in the canonical basis)
  sqrtMat[0] = eigVal[0]*eigVec[0]*eigVec[0] + eigVal[1]*eigVec[2]*eigVec[2];
  sqrtMat[1] = eigVal[0]*eigVec[0]*eigVec[1] + eigVal[1]*eigVec[2]*eigVec[3];
  sqrtMat[2] = eigVal[0]*eigVec[1]*eigVec[1] + eigVal[1]*eigVec[3]*eigVec[3];

  //--- Compute 1/sqrt(Lbd_M1)
  eigVal[0] = 1. / eigVal[0];
  eigVal[1] = 1. / eigVal[1];

  //--- Compute M1^(-1/2)  (in the canonical basis)
  sqrtMatInv[0] = eigVal[0]*eigVec[0]*eigVec[0] + eigVal[1]*eigVec[2]*eigVec[2];
  sqrtMatInv[1] = eigVal[0]*eigVec[0]*eigVec[1] + eigVal[1]*eigVec[2]*eigVec[3];
  sqrtMatInv[2] = eigVal[0]*eigVec[1]*eigVec[1] + eigVal[1]*eigVec[3]*eigVec[3];


  //--- M2bar = M1^-1/2 * M2 * M1^-1/2
  M2bar[0] = sqrtMatInv[0]*sqrtMatInv[0]*met2[0] 
           + 2*sqrtMatInv[0]*sqrtMatInv[1]*met2[1]
           + sqrtMatInv[1]*sqrtMatInv[1]*met2[2];
  M2bar[1] = sqrtMatInv[0]*sqrtMatInv[1]*met2[0] 
           + sqrtMatInv[1]*sqrtMatInv[1]*met2[1]
           + sqrtMatInv[0]*sqrtMatInv[2]*met2[1]
           + sqrtMatInv[1]*sqrtMatInv[2]*met2[2];
  M2bar[2] = sqrtMatInv[1]*sqrtMatInv[1]*met2[0] 
           + 2*sqrtMatInv[1]*sqrtMatInv[2]*met2[1]
           + sqrtMatInv[2]*sqrtMatInv[2]*met2[2];


  //--- Spectral decomposition of M2bar 
  ierr = DMPlexMetricReduction2d_Internal(M2bar, eigVal, eigVec);CHKERRQ(ierr);

  //    => normalized eigen vector ~> transformation matrix P 
  
  //    => We keep max(eigVal, 1)
  eigVal[0] = fmax(fabs(eigVal[0]), 1);
  eigVal[1] = fmax(fabs(eigVal[1]), 1);

  //--- We come back to the initial basis thanks to P and M1^(1/2)
  //-   Multiply by P = transpose(eigVec)  (cf storage by line vs column)
  Mintbar[0] = eigVal[0]*eigVec[0]*eigVec[0] + eigVal[1]*eigVec[2]*eigVec[2];
  Mintbar[1] = eigVal[0]*eigVec[0]*eigVec[1] + eigVal[1]*eigVec[2]*eigVec[3];
  Mintbar[2] = eigVal[0]*eigVec[1]*eigVec[1] + eigVal[1]*eigVec[3]*eigVec[3];

  //-   Multiply by M1^1/2
  metNew[0] = sqrtMat[0]*sqrtMat[0]*Mintbar[0] 
            + 2*sqrtMat[0]*sqrtMat[1]*Mintbar[1]
            + sqrtMat[1]*sqrtMat[1]*Mintbar[2];
  metNew[1] = sqrtMat[0]*sqrtMat[1]*Mintbar[0] 
            + sqrtMat[1]*sqrtMat[1]*Mintbar[1]
            + sqrtMat[0]*sqrtMat[2]*Mintbar[1]
            + sqrtMat[1]*sqrtMat[2]*Mintbar[2];
  metNew[2] = sqrtMat[1]*sqrtMat[1]*Mintbar[0] 
            + 2*sqrtMat[1]*sqrtMat[2]*Mintbar[1]
            + sqrtMat[2]*sqrtMat[2]*Mintbar[2];

  PetscFunctionReturn(0);
}



// transpose(V) * Met * V, with Met given as triangular upper matrix
static inline double product(double * V, double * Met) {
  double len = V[0]*V[0]*Met[0] + 2.*V[0]*V[1]*Met[1] + V[1]*V[1]*Met[2]; 
  return len;
}



#undef __FUNCT__
#define __FUNCT__ "DMPlexMetricGradation2d_Internal"
PetscErrorCode DMPlexMetricGradation2d_Internal(DM dm, PetscReal * metric, PetscReal * x, PetscReal * y) {
  
  PetscReal         beta = 1.4;
  
  PetscBool         correction;
  PetscInt        * verTag, eStart, eEnd, numEdges, vStart, vEnd, numVertices;
  PetscInt          iteCor, e, v, iVer1, iVer2, iMet1, iMet2, i;
  PetscReal         grownMet1[3], grownMet2[3], met1[3], met2[3], metNew1[3], metNew2[3], v12[2], v21[2];
  PetscReal         ln_beta, lengthEdge1, lengthEdge2, eta2_12, eta2_21, diff;
  const PetscInt  * cone;
  
  PetscInt          ierr;
  
  
  PetscFunctionBegin;

  ierr = DMPlexGetDepthStratum(dm, 1, &eStart, &eEnd);CHKERRQ(ierr);
  numEdges = eEnd - eStart;
//  printf("DEBUG  numEdges: %d\n", numEdges);
  ierr = DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd);CHKERRQ(ierr);
  numVertices = vEnd - vStart;

//  for (iVer1 = 0; iVer1 < numVertices; ++iVer1) {
//    printf("DEBUG  metric[iVer: %d] : %f %f %f %f\n", iVer1, metric[4*iVer1+0], metric[4*iVer1+1], metric[4*iVer1+2], metric[4*iVer1+3]);
//  }
  
  ierr = PetscMalloc1(numVertices, &verTag);CHKERRQ(ierr);
  for (v = 0; v < vEnd-vStart; ++v) verTag[v] = 1;
  
  ln_beta = log(beta);
 
//  printf("DEBUG  log(beta): %f\n", ln_beta);

  correction = PETSC_TRUE;
  iteCor = 0;
  while (correction && iteCor < 500) {
    iteCor++;
//    printf("EBUG  iteCor: %d\n", iteCor);
    correction = PETSC_FALSE;    
    for (e = eStart; e < eEnd; ++e) {

      ierr = DMPlexGetCone(dm, e, &cone);CHKERRQ(ierr);
      iVer1 = cone[0]-vStart; 
      iVer2 = cone[1]-vStart; 
      iMet1 = 4*iVer1; 
      iMet2 = 4*iVer2;

//	  printf("DEBUG  iEdg: %d,  iVer1, iVer2: %d %d, iMet1, iMet2: %d %d,  ver: [%f %f], [%f %f]\n", e-eStart, iVer1, iVer2, iMet1, iMet2, x[iVer1], y[iVer1], x[iVer2], y[iVer2]);
      
      if (verTag[iVer1] < iteCor && verTag[iVer2] < iteCor) continue;

      v12[0] = x[iVer2]-x[iVer1];
      v12[1] = y[iVer2]-y[iVer1];
      v21[0] = -v12[0];
      v21[1] = -v12[1];

      met1[0] = metric[iMet1]; met1[1] = metric[iMet1+1]; met1[2] = metric[iMet1+3]; 
      met2[0] = metric[iMet2]; met2[1] = metric[iMet2+1]; met2[2] = metric[iMet2+3]; 

//	  printf("DEBUG          v12: %f %f  met1: %f %f %f  met2: %f %f %f\n", v12[0], v12[1], met1[0], met1[1], met1[2], met2[0], met2[1], met2[2] );

      lengthEdge1 = sqrt(product(v12, met1)); 
      lengthEdge2 = sqrt(product(v21, met2));
      eta2_12 = 1+lengthEdge1*ln_beta;
      eta2_21 = 1+lengthEdge2*ln_beta;
      eta2_12 *= eta2_12;
      eta2_21 *= eta2_21;
      eta2_12 = 1./eta2_12;
      eta2_21 = 1./eta2_21;
//	  printf("DEBUG          lenEdg1,2 : %f %f   eta2_12,21: %f %f\n", lengthEdge1, lengthEdge2, eta2_12, eta2_21);
      for (i=0;i<3;++i) {
        grownMet1[i] = eta2_12*met1[i];
        grownMet2[i] = eta2_21*met2[i];
      }
      
      metricIntersection(met1, grownMet2, metNew1); // intersect M1 and eta12*M2
      metricIntersection(met2, grownMet1, metNew2);

//	  printf("DEBUG          grownMet1,2: %f %f %f   %f %f %f\n", grownMet1[0], grownMet1[1], grownMet1[2], grownMet2[0], grownMet2[1], grownMet2[2]);
//      printf("DEBUG          metNew1,2: %f %f %f   %f %f %f\n", metNew1[0], metNew1[1], metNew1[2], metNew2[0], metNew2[1], metNew2[2]);

      // compute norm met-metnew for each edge end
      // if norm > 1e-3, tag = 1, correction = PETSC_TRUE
      diff = fabs(met1[0]-metNew1[0])+fabs(met1[1]-metNew1[1])+fabs(met1[2]-metNew1[2]);
      diff /= (fabs(met1[0])+fabs(met1[1])+fabs(met1[2]));
      if (diff > 1.e-3) {
        metric[iMet1]   = metNew1[0];
        metric[iMet1+1] = metNew1[1];
        metric[iMet1+2] = metNew1[1];
        metric[iMet1+3] = metNew1[2];
        verTag[iVer1] = iteCor+1;
        correction = PETSC_TRUE;
      }
      diff = fabs(met2[0]-metNew2[0])+fabs(met2[1]-metNew2[1])+fabs(met2[2]-metNew2[2]);
      diff /= (fabs(met2[0])+fabs(met2[1])+fabs(met2[2]));
      if (diff > 1.e-3) {
        metric[iMet2]   = metNew2[0];
        metric[iMet2+1] = metNew2[1];
        metric[iMet2+2] = metNew2[1];
        metric[iMet2+3] = metNew2[2];
        verTag[iVer2] = iteCor+1;
        correction = PETSC_TRUE;
      }


    }
  }  


  PetscFree(verTag);
  PetscFunctionReturn(0);
}

