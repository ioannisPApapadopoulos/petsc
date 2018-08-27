#include <petsc/private/dmpleximpl.h>   /*I      "petscdmplex.h"   I*/




#undef __FUNCT__
#define __FUNCT__ "DMPlexMetricReduction3d_Internal"
// metric given as upper triangular matrix mat
PetscErrorCode DMPlexMetricReduction3d_Internal(PetscReal * met, PetscReal * eigVal, PetscReal * eigVec) 
{
  
  int    i, k, cas;
  double nrm, mat[6];
  double a, b, c, d, ap, bp, cp, alpha, beta, lbd1, lbd2, lbd3, delta, eps;
  double px, ppx, x0, x1, xmin, pmin, p[3], tmp; 
  double w1[3], w2[3], w3[3], v1[3], v2[3], v3[3];
  double v1nrm, v2nrm, v3nrm, w1nrm, w2nrm, w3nrm, vecNrm, nrmInv;
  
  double us6 = 1./6;
  double us3 = 1./3;
  
  PetscFunctionBegin;

  for (i=0; i<6; ++i) 
    mat[i] = met[i];

  nrm = fabs(mat[0]);
  for (i=1; i<6; ++i)
    nrm = fmax(nrm,fabs(mat[i]));
  
  // check for null matrices
  if ( nrm < 1e-100 ) {
    eigVal[0] = 0; eigVal[1] = 0; eigVal[2] = 0;
    for (i=0; i<9; ++i) eigVec[i] = 0;
    eigVec[0] = 1; eigVec[4] = 1; eigVec[8] = 1;  
    PetscFunctionReturn(0);
  }
  
  // normalize the matrix
  nrmInv = 1. / nrm;
  for (i=0; i<6; ++i) mat[i] *= nrmInv;
  

  a = -mat[0]-mat[3]-mat[5];			// = -Trace(mat)
  b = mat[0]*mat[3] + mat[0]*mat[5] + mat[3]*mat[5]
    - mat[1]*mat[1] - mat[2]*mat[2] - mat[4]*mat[4];
  c = mat[0]*(mat[4]*mat[4]-mat[3]*mat[5])	// = -Det(mat)
    + mat[1]*(mat[1]*mat[5]-mat[2]*mat[4])
    + mat[2]*(mat[2]*mat[3]-mat[1]*mat[4]);
  
  cas = 0;
  
  // P'(X) = ap*X^2 + bp*X + cp
  ap = 3;
  bp = 2*a;
  cp = b;
  
  // First look for double or triple roots  
  delta = bp*bp-4*ap*cp;
  eps   = bp*bp*1e-10;
  
  if ( delta > eps ) {	 
    
    // P' has two different roots:
    // if a root of P' is a root of P, it is a double root for P
    delta  = sqrt(delta);
    eigVal[0] = (-bp+delta)*us6; // first root of P'
    px = (((eigVal[0]+a)*eigVal[0])+b)*eigVal[0] +c;	
    if ( fabs(px) < 1e-15 ) {
      eigVal[1] = eigVal[0];
      eigVal[2] = -a - 2*eigVal[0]; // sum of eig. val. = Trace  (and -a = Trace)
      cas = 2;
      goto endEigVal;
    }
    eigVal[0] = (-bp-delta)*us6; // other root of P' 
    px = (((eigVal[0]+a)*eigVal[0])+b)*eigVal[0] +c;
    if ( fabs(px) < 1e-15 ) {
      eigVal[1] = eigVal[0];
      eigVal[2] = -a - 2*eigVal[0];
      cas = 2;
      goto endEigVal;
    }
    // else P has 3 single roots, see later
  }
  else if ( fabs(delta) <= eps ) { // delta = 0 
    // P' has a double root => P has a triple root    
    eigVal[0] = eigVal[1] = eigVal[2] = -bp*us6;
    cas = 3;
    goto endEigVal;
  }
  else {
    exit(11);   // P' cannot have no real root
  }
  
  // If P has 3 single roots, find the middle one with a Newton algorithm
  // the inflection point (ie P"(x)=0) is used as a starting point
  // TODO is this faster than using Cardano's formulas ?
  
  x0  = -a*us3;               // root of P''
  ppx = b-a*a*us3;            // P'(x0) 
  px = (2*a*a*us6+b)*a*us3+c;	// P(x0)
  
  xmin = x0;
  pmin = fabs(px);
  
  for (i=0; i<100; ++i) {
    
    x1  = x0-px/ppx;
    px = (((x1+a)*x1)+b)*x1+c;
    if ( fabs(px) < 1e-18 ) {
      eigVal[1] = x1;
      break;
    }
    if ( fabs(px) < pmin ) {
      xmin = x1;
      pmin = fabs(px);
    }
    ppx = (ap*x1+bp)*x1+cp; 
    if ( fabs((x1-x0)/x0) < 1e-7 ) {  // the algorithm has (almost) converged
      eigVal[1] = xmin;
      if ( pmin > 1e-12 ){ // check that P(xmin) really close to 0
        exit(22);
      }
      break;
    } 
    x0 = x1;
  }
  
  if ( i == 100 ) {
    eigVal[1] = xmin;
    if ( pmin > 1e-12 )
      exit(13);
  }
  
  
  // We have now P(X) = (X-eigVal[1])*(X^2 + alpha*X + beta)
  // so we just have to find the roots of quadratic polynomial
  alpha = a + eigVal[1];
  beta = b + eigVal[1]*alpha;
  
  delta = alpha*alpha-4*beta;
  eps   = alpha*alpha*1e-10; 
  
  if ( fabs(delta) < eps ) { // double root (lamda[0]=eigVal[1])
    eigVal[2] = eigVal[1];
    eigVal[0] = eigVal[1] = -0.5*alpha;
    cas = 2;
    goto endEigVal;
  }
  else if ( delta < 0. ) {
    exit(14);
  }
  
  delta  = sqrt(delta);
  eigVal[2] =  0.5*(delta-alpha); 
  eigVal[0] = -0.5*(delta+alpha); 
  
  //--- check another time if a double/triple roots is obtained at 1e-5
  //--- very important for the conditionning when looking for the eigenvectors
  for (i=0; i<3; ++i)
    p[i] = fabs( (((eigVal[i]+a)*eigVal[i] )+b)*eigVal[i]+c);
  lbd1 = fabs(eigVal[0]);
  lbd2 = fabs(eigVal[1]);
  lbd3 = fabs(eigVal[2]);
  
  if ( fabs(eigVal[0]-eigVal[1]) < 1e-5*(lbd1+lbd2)*0.5 ) {
    if ( p[0] < p[1] ) {
      eigVal[1] = eigVal[0];
      p[1] = p[0];
    }
    else {
      eigVal[0] = eigVal[1];
      p[0] = p[1];
    }
      
    // check for triple root
    if ( fabs(eigVal[1]-eigVal[2]) < 1e-5*(lbd2+lbd3)*0.5 ) {
      if ( p[1] < p[2] ) eigVal[2] = eigVal[1];
      else eigVal[1] = eigVal[2];
      cas = 3;
      goto endEigVal;
    }
    else 
      cas = 2;
      goto endEigVal;
  }
  
  else if ( fabs(eigVal[1]-eigVal[2]) < 1e-5*(lbd2+lbd3)*0.5 ) {
    if ( p[2] < p[1] ) {
      eigVal[1] = eigVal[2];
      p[1] = p[2];
    }
    else {
      eigVal[2] = eigVal[1];
      p[2] = p[1];
    }
    tmp = eigVal[0];	
    eigVal[0] = eigVal[2];
    eigVal[2] = tmp;    
    cas = 2;
    goto endEigVal;
  }
  
  else if ( fabs(eigVal[0]-eigVal[2]) < 1e-5*(lbd1+lbd3)*0.5 ) {
    if ( p[2] < p[0] ) {
      eigVal[0] = eigVal[2];
      p[0] = p[2];
    }
    else {
      eigVal[2] = eigVal[0];
      p[2] = p[0];
    }
    tmp = eigVal[1];
    eigVal[1] = eigVal[2];
    eigVal[2] = tmp;    
    cas = 2;
    goto endEigVal;
  }
  
  cas = 1;
  endEigVal:
  
  switch (cas) {
  
  // case with hree single eigenValues
  case 1:
    // vk= wi/\wj !=0  with  wi,wj line of W = mat-eigVal[i]*Id
    for (k=0; k<2; ++k) {
      w1[0] = mat[0]-eigVal[k]; w1[1] = mat[1];           w1[2] = mat[2];
      w2[0] = mat[1];           w2[1] = mat[3]-eigVal[k]; w2[2] = mat[4];
      w3[0] = mat[2];           w3[1] = mat[4];           w3[2] = mat[5]-eigVal[k];

      // cross products
      v1[0] =  w1[1]*w3[2] - w1[2]*w3[1]; v1[1] = -w1[0]*w3[2] + w1[2]*w3[0]; v1[2] =  w1[0]*w3[1] - w1[1]*w3[0];
      v1nrm = v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2];
      v2[0] =  w1[1]*w2[2] - w1[2]*w2[1]; v2[1] = -w1[0]*w2[2] + w1[2]*w2[0]; v2[2] =  w1[0]*w2[1] - w1[1]*w2[0];
      v2nrm = v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2];
      v3[0] =  w2[1]*w3[2] - w2[2]*w3[1]; v3[1] = -w2[0]*w3[2] + w2[2]*w3[0]; v3[2] =  w2[0]*w3[1] - w2[1]*w3[0];
      v3nrm = v3[0]*v3[0]+v3[1]*v3[1]+v3[2]*v3[2];

      // take the vector with the highest norm
      if ( v1nrm >= v2nrm && v1nrm >= v3nrm) {
       vecNrm = 1./sqrt(v1nrm);
        eigVec[3*k] = v1[0]*vecNrm; eigVec[3*k+1] = v1[1]*vecNrm; eigVec[3*k+2] = v1[2]*vecNrm;
      }
      else if ( v2nrm >= v1nrm && v2nrm >= v3nrm) {
       vecNrm = 1./sqrt(v2nrm);
        eigVec[3*k] = v2[0]*vecNrm; eigVec[3*k+1] = v2[1]*vecNrm; eigVec[3*k+2] = v2[2]*vecNrm;
      }
      else {
       vecNrm = 1./sqrt(v3nrm);
        eigVec[3*k] = v3[0]*vecNrm; eigVec[3*k+1] = v3[1]*vecNrm; eigVec[3*k+2] = v3[2]*vecNrm;
      }
    }

    // The last eigenvector is simply orthogonal to both others : v3=v1/\v2
    eigVec[6]  = eigVec[1]*eigVec[5] - eigVec[2]*eigVec[4];
    eigVec[7]  = eigVec[2]*eigVec[3] - eigVec[0]*eigVec[5];
    eigVec[8]  = eigVec[0]*eigVec[4] - eigVec[1]*eigVec[3];
   vecNrm = eigVec[6]*eigVec[6]+eigVec[7]*eigVec[7]+eigVec[8]*eigVec[8];
   vecNrm = 1./sqrt(vecNrm);
    eigVec[6] *=vecNrm; 
    eigVec[7] *=vecNrm; 
    eigVec[8] *=vecNrm; 
    
    break;

  // case with a double eigVal (eigVal[0] = eigVal[1]) and a single one (eigVal[2])  
  case 2:
    //  v3
    w1[0] = mat[0]-eigVal[2]; w1[1] = mat[1];           w1[2] = mat[2];
    w2[0] = mat[1];           w2[1] = mat[3]-eigVal[2]; w2[2] = mat[4];
    w3[0] = mat[2];           w3[1] = mat[4];           w3[2] = mat[5]-eigVal[2];

    // cross product of the lines
    v1[0] =  w1[1]*w3[2] - w1[2]*w3[1]; v1[1] = -w1[0]*w3[2] + w1[2]*w3[0]; v1[2] =  w1[0]*w3[1] - w1[1]*w3[0];
    v1nrm = v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2];
    v2[0] =  w1[1]*w2[2] - w1[2]*w2[1]; v2[1] = -w1[0]*w2[2] + w1[2]*w2[0]; v2[2] =  w1[0]*w2[1] - w1[1]*w2[0];
    v2nrm = v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2];
    v3[0] =  w2[1]*w3[2] - w2[2]*w3[1]; v3[1] = -w2[0]*w3[2] + w2[2]*w3[0]; v3[2] =  w2[0]*w3[1] - w2[1]*w3[0];
    v3nrm = v3[0]*v3[0]+v3[1]*v3[1]+v3[2]*v3[2];

    // take the vector with the highest norm
    if ( v1nrm >= v2nrm && v1nrm >= v3nrm) {
     vecNrm = 1./sqrt(v1nrm);
      eigVec[6] = v1[0]*vecNrm; eigVec[7] = v1[1]*vecNrm; eigVec[8] = v1[2]*vecNrm;
    }
    else if ( v2nrm >= v1nrm && v2nrm >= v3nrm) {
     vecNrm = 1./sqrt(v2nrm);
      eigVec[6] = v2[0]*vecNrm; eigVec[7] = v2[1]*vecNrm; eigVec[8] = v2[2]*vecNrm;
    }
    else {
     vecNrm = 1./sqrt(v3nrm);
      eigVec[6] = v3[0]*vecNrm; eigVec[7] = v3[1]*vecNrm; eigVec[8] = v3[2]*vecNrm;
    }    
    
    // v1
    w1nrm = w1[0]*w1[0]+w1[1]*w1[1]+w1[2]*w1[2];
    w2nrm = w2[0]*w2[0]+w2[1]*w2[1]+w2[2]*w2[2];
    w3nrm = w3[0]*w3[0]+w3[1]*w3[1]+w3[2]*w3[2];
    if ( w1nrm >=  w2nrm && w1nrm > w3nrm ) {
     vecNrm = 1./sqrt(w1nrm);
      eigVec[0] = w1[0]*vecNrm; eigVec[1] = w1[1]*vecNrm; eigVec[2] = w1[2]*vecNrm;
    }
    else if ( w2nrm >=  w1nrm && w2nrm > w3nrm ) {
     vecNrm = 1./sqrt(w2nrm);
      eigVec[0] = w2[0]*vecNrm; eigVec[1] = w2[1]*vecNrm; eigVec[2] = w2[2]*vecNrm;
    }
    else {
     vecNrm = 1./sqrt(w3nrm);
      eigVec[0] = w3[0]*vecNrm; eigVec[1] = w3[1]*vecNrm; eigVec[2] = w3[2]*vecNrm;
    }
        
    // The last eigenvector is simply orthogonal to both others: v2=v3/\v1
    eigVec[3]  = eigVec[7]*eigVec[2] - eigVec[8]*eigVec[1];
    eigVec[4]  = eigVec[8]*eigVec[0] - eigVec[6]*eigVec[2];
    eigVec[5]  = eigVec[6]*eigVec[1] - eigVec[7]*eigVec[0];
    vecNrm = eigVec[3]*eigVec[3]+eigVec[4]*eigVec[4]+eigVec[5]*eigVec[5];
    vecNrm = 1./sqrt(vecNrm);
    eigVec[3] *=vecNrm; eigVec[4] *=vecNrm; eigVec[5] *=vecNrm; 
    break;

  // triple eigenvalue => isotropic metric = lamba*Id  
  case 3:
    for (i=0; i<9; ++i) 
      eigVec[i] = 0;
    eigVec[0] = 1; eigVec[4] = 1; eigVec[8] = 1;  
    break;  
  
  default:
    exit(15);
  }
  
  eigVal[0] *= nrm;
  eigVal[1] *= nrm;
  eigVal[2] *= nrm;
  
  PetscFunctionReturn(0);
  
}







static inline void getMetric3(double* eigVal, double* eigVec, double* met)
{
  met[0] = eigVal[0]*eigVec[0]*eigVec[0] + eigVal[1]*eigVec[3]*eigVec[3] + eigVal[2]*eigVec[6]*eigVec[6];
  met[1] = eigVal[0]*eigVec[0]*eigVec[1] + eigVal[1]*eigVec[3]*eigVec[4] + eigVal[2]*eigVec[6]*eigVec[7];
  met[2] = eigVal[0]*eigVec[0]*eigVec[2] + eigVal[1]*eigVec[3]*eigVec[5] + eigVal[2]*eigVec[6]*eigVec[8];
  met[3] = eigVal[0]*eigVec[1]*eigVec[1] + eigVal[1]*eigVec[4]*eigVec[4] + eigVal[2]*eigVec[7]*eigVec[7];
  met[4] = eigVal[0]*eigVec[1]*eigVec[2] + eigVal[1]*eigVec[4]*eigVec[5] + eigVal[2]*eigVec[7]*eigVec[8];
  met[5] = eigVal[0]*eigVec[2]*eigVec[2] + eigVal[1]*eigVec[5]*eigVec[5] + eigVal[2]*eigVec[8]*eigVec[8];
}

/*
  Matrix product of 2 matrix of type "transformation matrix": P'LP or PLP'
  P is ordinary, is L is diagonal, the result is symmetric
      [ p0 p1 p2 ]        [L0 0 0]          [ res0 res1 res2 ] 
  P = [ p3 p4 p5 ]    L = [0 L1 0]    res = [  -   res3 res4 ] 
      [ p6 p7 p8 ]        [0 0 L2]          [  -    -   res5 ] 
*/

static inline void matrixProduct_PLtP(double* P, double* L, double* res)
{
  res[0] = L[0]*P[0]*P[0] + L[1]*P[1]*P[1] + L[2]*P[2]*P[2];
  res[1] = L[0]*P[0]*P[3] + L[1]*P[1]*P[4] + L[2]*P[2]*P[5];
  res[2] = L[0]*P[0]*P[6] + L[1]*P[1]*P[7] + L[2]*P[2]*P[8];
  res[3] = L[0]*P[3]*P[3] + L[1]*P[4]*P[4] + L[2]*P[5]*P[5];
  res[4] = L[0]*P[3]*P[6] + L[1]*P[4]*P[7] + L[2]*P[5]*P[8];
  res[5] = L[0]*P[6]*P[6] + L[1]*P[7]*P[7] + L[2]*P[8]*P[8]; 
}



/*
  Matrix product of 2 matrix of type "transformation matrix": SMS
  S, M and the result are symmetric
      [ s0 s1 s2 ]        [m0 m1 m2]          [ res0 res1 res2 ] 
  S = [ -  s3 s4 ]    M = [ - m3 m4]    res = [  -   res3 res4 ] 
      [ -  -  s5 ]        [-  -  m5]          [  -    -   res5 ]
*/
static inline void matrixProduct_SMS(double* S, double* M, double* res)
{  
  double mat[9];
  
  //--- First mat=SM
  mat[0] = S[0]*M[0] + S[1]*M[1] + S[2]*M[2];
  mat[1] = S[0]*M[1] + S[1]*M[3] + S[2]*M[4];
  mat[2] = S[0]*M[2] + S[1]*M[4] + S[2]*M[5];
  mat[3] = S[1]*M[0] + S[3]*M[1] + S[4]*M[2];
  mat[4] = S[1]*M[1] + S[3]*M[3] + S[4]*M[4];
  mat[5] = S[1]*M[2] + S[3]*M[4] + S[4]*M[5];
  mat[6] = S[2]*M[0] + S[4]*M[1] + S[5]*M[2];
  mat[7] = S[2]*M[1] + S[4]*M[3] + S[5]*M[4];
  mat[8] = S[2]*M[2] + S[4]*M[4] + S[5]*M[5];
  

  //--- Second res=SMS
  res[0] = mat[0]*S[0] + mat[1]*S[1] + mat[2]*S[2];
  res[1] = mat[0]*S[1] + mat[1]*S[3] + mat[2]*S[4];
  res[2] = mat[0]*S[2] + mat[1]*S[4] + mat[2]*S[5];
  res[3] = mat[3]*S[1] + mat[4]*S[3] + mat[5]*S[4];
  res[4] = mat[3]*S[2] + mat[4]*S[4] + mat[5]*S[5];
  res[5] = mat[6]*S[2] + mat[7]*S[4] + mat[8]*S[5];
   
}

/*

  Application of a metric to a vector: res = <Mv, v> = v'Mv
  v = [v0, v1, v2]'          [m0 m1 m2]
  M is symmetric matrix M =  [-  m3 m4]
                             [-  -  m5]
  res is a scalar
  
*/
static inline void matrixProduct_tvMv(double* v, double* M, double* res)
{
  *res = v[0]*(v[0]*M[0] + v[1]*M[1] + v[2]*M[2]) +
         v[1]*(v[0]*M[1] + v[1]*M[3] + v[2]*M[4]) +
         v[2]*(v[0]*M[2] + v[1]*M[4] + v[2]*M[5]) ;
}


static inline double getMatDet3(double* mat)
{
  double aa,bb,cc,det;
  
  aa = mat[4]*mat[8] - mat[5]*mat[7];
  bb = mat[5]*mat[6] - mat[3]*mat[8];
  cc = mat[3]*mat[7] - mat[4]*mat[6];
  
  det = mat[0]*aa + mat[1]*bb + mat[2]*cc; 

  return det;
}

static inline double getMatInv3(double *mat, double *matInv)
{
  double det,detInv;
  
  det = getMatDet3(mat);
  if ( fabs(det) < 1e-100 ) {
    fprintf(stderr,"\n  ## ERROR getMatInv : An non-inversible matrix is treated !\n");
    fprintf(stderr,"	        [%12.8e  %12.8e  %12.8e]\n",mat[0],mat[1],mat[2]);
    fprintf(stderr,"     M  = [%12.8e  %12.8e  %12.8e]  with  det = %12.8e\n",mat[3],mat[4],mat[5],det);
    fprintf(stderr,"          [%12.8e  %12.8e  %12.8e]\n",mat[6],mat[7],mat[8]);
    exit(1);
  }
  detInv = 1./det;
  
  matInv[0] = (mat[4]*mat[8] - mat[5]*mat[7])*detInv;
  matInv[1] = (mat[2]*mat[7] - mat[1]*mat[8])*detInv;
  matInv[2] = (mat[1]*mat[5] - mat[2]*mat[4])*detInv;
  matInv[3] = (mat[5]*mat[6] - mat[3]*mat[8])*detInv;
  matInv[4] = (mat[0]*mat[8] - mat[2]*mat[6])*detInv;
  matInv[5] = (mat[2]*mat[3] - mat[0]*mat[5])*detInv;
  matInv[6] = (mat[3]*mat[7] - mat[4]*mat[6])*detInv;
  matInv[7] = (mat[1]*mat[6] - mat[0]*mat[7])*detInv;
  matInv[8] = (mat[0]*mat[4] - mat[1]*mat[3])*detInv;

  return det;
}



static inline void sortEigen(double *eigVal, double *eigVec)
{
  int    i;  
  double eigValTmp[3],eigVecTmp[9];

  //--- Init
  for (i=0; i<3; ++i)
    eigValTmp[i] = eigVal[i];
  for (i=0; i<9; ++i)
    eigVecTmp[i] = eigVec[i];

  //--- Switch
  if ( eigVal[0] > eigVal[1] ) {        // 0 > 1
    if ( eigVal[1] > eigVal[2] ) {      // 0 > 1 > 2
      // already sorted
    }
    else if ( eigVal[2] > eigVal[0] ) { // 2 > 0 > 1
      eigVal[0] = eigValTmp[2];
      eigVal[1] = eigValTmp[0]; 
      eigVal[2] = eigValTmp[1]; 
      for (i=0; i<3; ++i) {
        eigVec[0+i] = eigVecTmp[6+i];
        eigVec[3+i] = eigVecTmp[0+i];
        eigVec[6+i] = eigVecTmp[3+i];
      }       
    }
    else {                              // 0 > 2 > 1
      eigVal[1] = eigValTmp[2];
      eigVal[2] = eigValTmp[1];
      for (i=0; i<3; ++i) {
        eigVec[3+i] = eigVecTmp[6+i];
        eigVec[6+i] = eigVecTmp[3+i];
      }       
    } 
  }
  
  else {                                // 1 > 0
    if ( eigVal[0] > eigVal[2] ) {      // 1 > 0 > 2
      eigVal[0] = eigValTmp[1];
      eigVal[1] = eigValTmp[0]; 
      for (i=0; i<3; ++i) {
        eigVec[0+i] = eigVecTmp[3+i];
        eigVec[3+i] = eigVecTmp[0+i];
      }       
    }
    else if ( eigVal[2] > eigVal[1] ) { // 2 > 1 > 0
      eigVal[0] = eigValTmp[2];
      eigVal[2] = eigValTmp[0]; 
      for (i=0; i<3; ++i) {
        eigVec[0+i] = eigVecTmp[6+i];
        eigVec[6+i] = eigVecTmp[0+i];
      }       
    }
    else {                              // 1 > 2 > 0
      eigVal[0] = eigValTmp[1];
      eigVal[1] = eigValTmp[2];
      eigVal[2] = eigValTmp[0];        
      for (i=0; i<3; ++i) {
        eigVec[0+i] = eigVecTmp[3+i];
        eigVec[3+i] = eigVecTmp[6+i];
        eigVec[6+i] = eigVecTmp[0+i];
      }       
    } 
  }
  
}


/*

  Case (e11 e12 0), (e21 e22 0) : intersect 2 infinite elliptic-cylinders
  Case 4. in the report
  
  There are 2 subcases: a. The degenerate directions coincide
                        b. Otherwise
                        
  Assume the eigen values are ordered from the largest to the lowest
                        
*/
static inline void degenerateIntersection_11(double * met1, double *eigVec1, 
                                                  double * met2, double *eigVec2, double * metNew)
{  
  int    i;
  double vec[3],nrm,Pmat[9],Pinv[9];
  double eigVal[3],siz1,siz2;

  double coef = sqrt(3)/2.;
  
  //--- Cross product of the two degenerate directions: w1/\w2
  vec[0] = eigVec1[7]*eigVec2[8] - eigVec1[8]*eigVec2[7];
  vec[1] = eigVec1[8]*eigVec2[6] - eigVec1[6]*eigVec2[8];
  vec[2] = eigVec1[6]*eigVec2[7] - eigVec1[7]*eigVec2[6];
  nrm    = vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2];
  
  //--- Check subcases
  if ( nrm < 1e-30 ) {
    //- case 4.(a)  in the report
    fprintf(stderr, "  Metric intersection: degenerate case not treated yet: case 4.(a)\n");
    //- hack for robustness: my accurate approx
    for (i=0; i<6; ++i)
      metNew[i] = coef*(met1[i]+met2[i]);
  }
  
  else {
    //--- case 4.(b)  in the report: compute the metric in basis (w1,w2,w1/\w2=vec)
    //- Compute w1="M1 degenerate direction" length in M2 : <w1,M2w1>
    matrixProduct_tvMv(&eigVec1[6], met2, &eigVal[0]);

    //- Compute w2="M2 degenerate direction" length in M1 : <w2,M1w2>
    matrixProduct_tvMv(&eigVec2[6], met1, &eigVal[1]);

    //- Perform the 1D metric intersection in the orthogonal direction: max(vec*M1*vec,vec*M2*vec)
    matrixProduct_tvMv(vec, met1, &siz1);
    matrixProduct_tvMv(vec, met2, &siz2);
    eigVal[2] = fmax(siz1,siz2);
  
    //- Set transformation matrix P and compute M = tP^-1 * L * P^-1 
    Pmat[0] = eigVec1[6]; Pmat[1] = eigVec1[7]; Pmat[2] = eigVec1[8];
    Pmat[3] = eigVec2[6]; Pmat[4] = eigVec2[7]; Pmat[5] = eigVec2[8];
    Pmat[6] = vec[0];     Pmat[7] = vec[1];     Pmat[8] = vec[2];
    getMatInv3(Pmat, Pinv);
    
    matrixProduct_PLtP(Pinv, eigVal, metNew);
  }
  
}







/*

  Case (e11 0 0), (e21 0 0) : intersect 2 infinite board-shapes
  Case 2. in the report
  
  There are 2 subcases: a. The NON-degenerate directions coincide ie parallel plane
                        b. Otherwise
                        
  Assume the eigen values are ordered from the largest to the lowest
                        
*/
static inline void degenerateIntersection_22(double * met1, double * eigVal1, double *eigVec1, 
                                              double * met2, double * eigVal2, double *eigVec2, double * metNew)
{  
  int    i;
  double vec[3],nrm;
  double eigVal[3],siz1,siz2;

  double coef = sqrt(3)/2.;


  //--- Cross product of the two NON-degenerate directions: u1/\u2
  vec[0] = eigVec1[1]*eigVec2[2] - eigVec1[2]*eigVec2[1];
  vec[1] = eigVec1[2]*eigVec2[0] - eigVec1[0]*eigVec2[2];
  vec[2] = eigVec1[0]*eigVec2[1] - eigVec1[1]*eigVec2[0];
  nrm    = vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2];

  //--- Check subcases
  if ( nrm < 1.e-30 ) {
    //- case 2.(a)  in the report: compute the metric in basis (u1,v1,w1)
    //- Perform the 1D metric intersection in the NON-degenerate direction: max(vec*M1*vec,vec*M2*vec)=max(l1,l2)
    siz1 = eigVal1[0];  //Spy_MatrixProduct_tvMv(vec, met1, &siz1);
    siz2 = eigVal2[0];  //Spy_MatrixProduct_tvMv(vec, met2, &siz2);
    eigVal[0] = fmax(siz1,siz2);

    //- The other 2 directions are degenerate => eigVal=0
    eigVal[1] = eigVal[2] = 0.;
    
    //- Take one of the possible basis, here eigVec1
    getMetric3(eigVal,eigVec1,metNew);
  }
  
  else {
    //--- case 2.(b)  in the report
    fprintf(stderr, "  Metric intersection: degenerate case not treated yet: case 2.(b)\n");
    //- hack for robustness: my accurate approx
    for (i=0; i<6; ++i)
      metNew[i] = coef*(met1[i]+met2[i]);

    //TODO: Pas trop dur
  }
}





/*

  Case (e11 0 0), (e21 e22 0) : intersect an infinite elliptic-cylinder and an infinite board-shape
  Case 3. in the report
  
  There are 2 subcases: a. The non-degenerate direction u1 coincides with degenerate direction w2
                        b. Otherwise
                        
  Assume the eigen values are ordered from the largest to the lowest
                        
*/
static inline void degenerateIntersection_21(double * met1, double *eigVec1, 
                                              double * met2, double *eigVec2, double * metNew)
{  
  int    i;
  double sca,Pmat[9],Pinv[9];
  double eigVal[3];

  double coef = sqrt(3)/2.;


  //--- Dot product of the non-degenerate direction and the degenerate one: u1.w2
  sca = eigVec1[0]*eigVec2[6] + eigVec1[1]*eigVec2[7] + eigVec1[2]*eigVec2[8];

  //--- Check subcases
  if ( sca < 1e-30 ) {
    //- case 3.(a)  in the report
    fprintf(stderr, "  Metric intersection: degenerate case not treated yet: case 3.(a)\n");
    //- hack for robustness: my accurate approx
    for (i=0; i<6; ++i)
      metNew[i] = coef*(met1[i]+met2[i]);
  }
  
  else {
    //--- case 3.(b)  in the report: compute the metric in degenerate direction basis (v1,w1,w2)
    //- Compute v1/w1="M1 degenerate direction" length in M2 : <w1,M2w1>
    matrixProduct_tvMv(&eigVec1[3], met2, &eigVal[0]);
    matrixProduct_tvMv(&eigVec1[6], met2, &eigVal[1]);

    //- Compute w2="M2 degenerate direction" length in M1 : <w2,M1w2>
    matrixProduct_tvMv(&eigVec2[6], met1, &eigVal[2]);

    //- Set transformation matrix P and compute M = tP^-1 * L * P^-1 
    Pmat[0] = eigVec1[3]; Pmat[1] = eigVec1[4]; Pmat[2] = eigVec1[5];
    Pmat[3] = eigVec1[6]; Pmat[4] = eigVec1[7]; Pmat[5] = eigVec1[8];
    Pmat[6] = eigVec2[6]; Pmat[7] = eigVec2[7]; Pmat[8] = eigVec2[8];
    getMatInv3(Pmat, Pinv);
    
    matrixProduct_PLtP(Pinv, eigVal, metNew);
  }
}







static inline void intersectdegenerateCase3D(double *met1, double *eigVal1, double *eigVec1, 
                                           double *met2, double *eigVal2, double *eigVec2,
                                           int    *cpt,  double *metNew)
{
  
  
  //--- Sort eigenValues from largest to lowest
  sortEigen(eigVal1,eigVec1);
  sortEigen(eigVal2,eigVec2);
   
  //--- Dissociated cases depending on the number of null eigenvalues
  if ( cpt[0] == 1 ) {
    if ( cpt[1] == 1 ) {
      //--- case (e11 e12 0), (e21 e22 0) : intersect 2 infinite elliptic-cylinders
      degenerateIntersection_11(met1,eigVec1,met2,eigVec2,metNew);
    }
    else {
      //--- case (e11 e12 0), (e21, 0 0) : intersect an infinite elliptic-cylinder and an infinite board-shape
      degenerateIntersection_21(met2,eigVec2,met1,eigVec1,metNew);
    }
  }
  else {
    if ( cpt[1] == 1 ) {
      //--- case (e11 0 0), (e21 e22 0) : intersect an infinite elliptic-cylinder and an infinite board-shape
      degenerateIntersection_21(met1,eigVec1,met2,eigVec2,metNew);
    }
    else {
      //--- case (e11 0 0), (e21, 0 0) : intersect 2 infinite board-shapes
      degenerateIntersection_22(met1,eigVal1,eigVec1,met2,eigVal2,eigVec2,metNew);
    }
  }
}


PetscErrorCode metricIntersection3d(double* metric1, double* metric2, double* metNew)
{
  
  int    i,cpt[2];
  double det1,det2,tmp;
  double met1[6],met2[6],eigVal[3],eigVec[9], eigVal1[3],eigVec1[9], eigVal2[3],eigVec2[9];
  double sqrtM1[6], sqrtInvM1[6], M2bar[6], metNewBar[6];
  
  PetscInt ierr;

  double detMin = 1e-10;
  double eigMin = 1e-15;
  
  PetscFunctionBegin;
  
  //--- Copy metric in working array because they may be switched
  for (i=0; i<6; ++i) {    
    met1[i] = metric1[i]; 
    met2[i] = metric2[i]; 
  }

//  printf("DEBUG      met1: %f %f %f   met2: %f %f %f\n", met1[0], met1[3], met1[5], met2[0], met2[3], met2[5]);
  

  //--- First check for null metric
  det1 = met1[0] * ( met1[3]*met1[5] - met1[4]*met1[4]) - met1[1] * ( met1[1]*met1[5] - met1[2]*met1[4]) + met1[2] * ( met1[1]*met1[4] - met1[2]*met1[3]);;
  det2 = met2[0] * ( met2[3]*met2[5] - met2[4]*met2[4]) - met2[1] * ( met2[1]*met2[5] - met2[2]*met2[4]) + met2[2] * ( met2[1]*met2[4] - met2[2]*met2[3]);
  if ( det1 < detMin ) {
    if ( det2 < detMin ){
      //--- Both metric are zero => metNew = 0
      for (i=0; i<6; ++i) metNew[i] = 0.;
      return 3;      
    }
    else {
      for (i=0; i<6; ++i) metNew[i] = met2[i];
      return 3;
    }
  }
  if ( det2 < detMin ) {
    for (i=0; i<6; ++i) metNew[i] = met1[i];
    return 3;
  }
  

  //--- Spectral decompistion of the two metric
  ierr = DMPlexMetricReduction3d_Internal(met1,eigVal1,eigVec1);CHKERRQ(ierr);  
  ierr = DMPlexMetricReduction3d_Internal(met2,eigVal2,eigVec2);CHKERRQ(ierr);
  
//  printf("DEBUG        eigVal1: %f %f %f, eigVec1: %f %f %f  %f %f %f  %f %f %f\n", eigVal1[0], eigVal1[1], eigVal1[2],
//    eigVec1[0], eigVec1[1], eigVec1[2], eigVec1[3], eigVec1[4], eigVec1[5], eigVec1[6], eigVec1[7], eigVec1[8]);
//  printf("DEBUG        eigVal2: %f %f %f, eigVec2: %f %f %f  %f %f %f  %f %f %f\n", eigVal2[0], eigVal2[1], eigVal2[2],
//    eigVec2[0], eigVec2[1], eigVec2[2], eigVec2[3], eigVec2[4], eigVec2[5], eigVec2[6], eigVec2[7], eigVec2[8]);

  
  //--- check for degenerate cases
  // -0- if M1 and M2 def > 0   => default
  // -1- M1 and M2 are zeros    => metNew = 0
  // -2- M1 == 0                => metNew = M2
  // -3- M2 == 0                => metNew = M1
  // -4- M1 def > 0 and M2 not  => default
  // -5- M2 def > 0 and M1 not  => switch M1 and M2
  // -6- M1 and M2 not def > 0  => welcome in hell
  
  cpt[0] = cpt[1] = 0;
  for (i=0; i<3; ++i) {
    if ( fabs(eigVal1[i]) < eigMin ) cpt[0]++;
    if ( fabs(eigVal2[i]) < eigMin ) cpt[1]++;
  }
  
  if ( cpt[0] == 0 ) {
    if ( cpt[1] != 3 ) {    // case 0 and 4
      //- do nothing, normal case
      for (i=0; i<3; ++i) eigVal[i] = eigVal1[i];
      for (i=0; i<9; ++i) eigVec[i] = eigVec1[i];
    }
    else {                  // case 3
      for (i=0; i<6; ++i) metNew[i] = met1[i];
      return 3;           
    }
  }
  else if ( cpt[0] == 1 || cpt[0] == 2 ) {
    if ( cpt[1] == 0 ) {
      //- switch M1 and M2
      for (i=0; i<6; ++i){
        tmp     = met1[i];
        met1[i] = met2[i];
        met2[i] = tmp;
      }
      for (i=0; i<3; ++i) eigVal[i] = eigVal2[i];
      for (i=0; i<9; ++i) eigVec[i] = eigVec2[i];
    }
    else if ( cpt[0] == 1 || cpt[0] == 2 ) {
      //--- Treat specific degenerate case
      intersectdegenerateCase3D(met1,eigVal1,eigVec1,met2,eigVal2,eigVec2,cpt,metNew);
      return 3;
    }
    else { // case 3
      for (i=0; i<6; ++i) metNew[i] = met1[i];
      return 3;           
    }
  }
  else {
    if ( cpt[1] != 3 ) { // case 2
      for (i=0; i<6; ++i) metNew[i] = met2[i];
      return 3;           
    }
    else {   // case 1
      for (i=0; i<6; ++i) metNew[i] = 0.;
      return 3;           
    } 
  }

  
  //--- Classic metric intersection with at least one non degenerate metric
  //--- Note that eigVec is the transpose of the transformation matrix P
  //--- (with P : canonical basis -> diagonalisation basis)
  //--- so sqrt(M1) = P sqrt(diag(eigVal1)) P' = eigVec1' sqrt(diag(eigVal1)) eigVec1
  
  //- Compute sqrt(Lbd_M1) ~> M1^(1/2)  (in the canonical basis)
  for (i=0; i<3;++i) 
    eigVal[i] = sqrt(fabs(eigVal[i]));
  getMetric3(eigVal,eigVec,sqrtM1);

//  printf("DEBUG        sqrtM1: %f %f %f %f %f %f\n", sqrtM1[0], sqrtM1[1], sqrtM1[2], sqrtM1[3], sqrtM1[4], sqrtM1[5]);
  

  //- Compute 1./sqrt(Lbd_M1) ~> M1^(-1/2)  (in the canonical basis)
  for (i=0; i<3;++i) 
    eigVal[i] = 1./eigVal[i];
  getMetric3(eigVal,eigVec,sqrtInvM1);

//  printf("DEBUG        sqrtInvM1: %f %f %f %f %f %f\n", sqrtInvM1[0], sqrtInvM1[1], sqrtInvM1[2], sqrtInvM1[3], sqrtInvM1[4], sqrtInvM1[5]);

  
  //--- Compute M2bar = M1^(-1/2) * M2 * M1^(-1/2)
  matrixProduct_SMS(sqrtInvM1,met2,M2bar);
//  printf("DEBUG        M2bar: %f %f %f %f %f %f\n", M2bar[0], M2bar[1], M2bar[2], M2bar[3], M2bar[4], M2bar[5]);

  //--- Spectral decomposition of M2bar 
  ierr = DMPlexMetricReduction3d_Internal(M2bar,eigVal,eigVec);CHKERRQ(ierr);


  //--- Intersection M1 and M2 comes to intersect M2bar with the unit circle 
  //--- and finally  M1/\M2 = M1^(1/2) * M2Bar * M1^(1/2)
  for (i=0; i<3; ++i)
    eigVal[i] = fmax(fabs(eigVal[i]),1.);
  getMetric3(eigVal,eigVec,metNewBar);

  matrixProduct_SMS(sqrtM1,metNewBar,metNew);

  return(0);
  
  
}



static inline double product3d(double * V, double * Met) {
  double len = V[0]*V[0]*Met[0] + V[1]*V[1]*Met[3] + V[2]*V[2]*Met[5] 
      + 2.*(V[0]*V[1]*Met[1] + V[0]*V[2]*Met[2] + V[1]*V[2]*Met[4]);
  return len;
}










#undef __FUNCT__
#define __FUNCT__ "DMPlexMetricGradation2d_Internal"
PetscErrorCode DMPlexMetricGradation3d_Internal(DM dm, PetscReal * metric, PetscReal * x, PetscReal * y, PetscReal * z) {
  
  PetscReal         beta = 1.4;
  
  PetscBool         correction;
  PetscInt        * verTag, eStart, eEnd, numEdges, vStart, vEnd, numVertices;
  PetscInt          iteCor, e, v, iVer1, iVer2, iMet1, iMet2, i;
  PetscReal         grownMet1[6], grownMet2[6], met1[6], met2[6], metNew1[6], metNew2[6], v12[3], v21[3];
  PetscReal         ln_beta, lengthEdge1, lengthEdge2, eta2_12, eta2_21, diff;
  const PetscInt  * cone;
  
  PetscInt          ierr;
  
  
  PetscFunctionBegin;

  ierr = DMPlexGetDepthStratum(dm, 1, &eStart, &eEnd);CHKERRQ(ierr);
  numEdges = eEnd - eStart;
  printf("DEBUG  numEdges: %d\n", numEdges);
  ierr = DMPlexGetDepthStratum(dm, 0, &vStart, &vEnd);CHKERRQ(ierr);
  numVertices = vEnd - vStart;
  
  ierr = PetscMalloc1(numVertices, &verTag);CHKERRQ(ierr);
  for (v = 0; v < vEnd-vStart; ++v) verTag[v] = 1;
  
  ln_beta = log(beta);
 
  correction = PETSC_TRUE;
  iteCor = 0;
  while (correction && iteCor < 500) {
//    printf("DEBUG  iteCor: %d\n", iteCor);
    iteCor++;
    correction = PETSC_FALSE;
    for (e = eStart; e < eEnd; ++e) {

      ierr = DMPlexGetCone(dm, e, &cone);CHKERRQ(ierr);
      iVer1 = cone[0]-vStart; 
      iVer2 = cone[1]-vStart; 
      iMet1 = 9*iVer1; 
      iMet2 = 9*iVer2;
      
//      printf("DEBUG  iEdg: %d,  iVer1, iVer2: %d %d,  ver: [%f %f %f], [%f %f %f]\n", e-eStart, iVer1, iVer2, x[iVer1], y[iVer1], z[iVer1], x[iVer2], y[iVer2], z[iVer2]);

      if (verTag[iVer1] < iteCor && verTag[iVer2] < iteCor) continue;

      v12[0] = x[iVer2]-x[iVer1];
      v12[1] = y[iVer2]-y[iVer1];
      v12[2] = z[iVer2]-z[iVer1];
      v21[0] = -v12[0];
      v21[1] = -v12[1];
      v21[2] = -v12[2];

      met1[0] = metric[iMet1]; met1[1] = metric[iMet1+1]; met1[2] = metric[iMet1+2];
                               met1[3] = metric[iMet1+4]; met1[4] = metric[iMet1+5];
                                                          met1[5] = metric[iMet1+8];
      met2[0] = metric[iMet2]; met2[1] = metric[iMet2+1]; met2[2] = metric[iMet2+2];
                               met2[3] = metric[iMet2+4]; met2[4] = metric[iMet2+5];
                                                          met2[5] = metric[iMet2+8];

//      printf("DEBUG      met1: %f %f %f   met2: %f %f %f\n", met1[0], met1[3], met1[5], met2[0], met2[3], met2[5]);

      lengthEdge1 = sqrt(product3d(v12, met1)); 
      lengthEdge2 = sqrt(product3d(v21, met2));
      eta2_12 = 1+lengthEdge1*ln_beta;
      eta2_21 = 1+lengthEdge2*ln_beta;
      eta2_12 *= eta2_12;
      eta2_21 *= eta2_21;
      eta2_12 = 1./eta2_12;
      eta2_21 = 1./eta2_21;
      for (i=0;i<6;++i) {
        grownMet1[i] = eta2_12*met1[i];
        grownMet2[i] = eta2_21*met2[i];
      }

//      printf("DEBUG      grownMet1: %f %f %f   grownMet2: %f %f %f\n", grownMet1[0], grownMet1[3], grownMet1[5], grownMet2[0], grownMet2[3], grownMet2[5]);
      
      metricIntersection3d(met1, grownMet2, metNew1); // intersect M1 and eta12*M2
      metricIntersection3d(met2, grownMet1, metNew2);

//      printf("DEBUG      metNew1: %f %f %f   metNew2: %f %f %f\n", metNew1[0], metNew1[3], metNew1[5], metNew2[0], metNew2[3], metNew2[5]);


      // compute norm met-metnew for each edge end
      // if norm > 1e-3, tag = 1, correction = PETSC_TRUE
      diff = 0;
      for (i=0; i<6; ++i) 
        diff += fabs(met1[i]-metNew1[i]);
      diff /= (fabs(met1[0])+fabs(met1[1])+fabs(met1[2])+fabs(met1[3])+fabs(met1[4])+fabs(met1[5]));
      if (diff > 1.e-3) {
        metric[iMet1]   = metNew1[0]; metric[iMet1+1] = metNew1[1]; metric[iMet1+2] = metNew1[2];
        metric[iMet1+3] = metNew1[1]; metric[iMet1+4] = metNew1[3]; metric[iMet1+5] = metNew1[4];
        metric[iMet1+6] = metNew1[2]; metric[iMet1+7] = metNew1[4]; metric[iMet1+8] = metNew1[5];
        verTag[iVer1] = iteCor+1;
        correction = PETSC_TRUE;
      }
      diff = 0;
      for (i=0; i<6; ++i) 
        diff += fabs(met2[i]-metNew2[i]);
      diff /= (fabs(met2[0])+fabs(met2[1])+fabs(met2[2])+fabs(met2[3])+fabs(met2[4])+fabs(met2[5]));
      if (diff > 1.e-3) {
        metric[iMet2]   = metNew2[0]; metric[iMet2+1] = metNew2[1]; metric[iMet2+2] = metNew2[2];
        metric[iMet2+3] = metNew2[1]; metric[iMet2+4] = metNew2[3]; metric[iMet2+5] = metNew2[4];
        metric[iMet2+6] = metNew2[2]; metric[iMet2+7] = metNew2[4]; metric[iMet2+8] = metNew2[5];
        verTag[iVer2] = iteCor+1;
        correction = PETSC_TRUE;
      }
      
//      if (e-eStart > 2) break;

    }
//    break;
  }  


  PetscFree(verTag);
  PetscFunctionReturn(0);
  
}
