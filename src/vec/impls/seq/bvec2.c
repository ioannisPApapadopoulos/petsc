

#ifndef lint
static char vcid[] = "$Id: bvec2.c,v 1.10 1995/03/17 04:55:42 bsmith Exp bsmith $";
#endif
/*
   Defines the sequential BLAS based vectors
*/

#include "sys/flog.h"
#include "inline/dot.h"
#include "inline/vmult.h"
#include "inline/setval.h"
#include "inline/copy.h"
#include "inline/axpy.h"
#include <math.h>
#include "vecimpl.h" 
#include "dvecimpl.h" 

#include "../bvec1.c"
#include "../dvec2.c"

static int VeiDVBCreateVector(Vec,Vec*);

static struct _VeOps DvOps = { VeiDVBCreateVector, 
              Veiobtain_vectors, Veirelease_vectors, VeiDVBdot,
              VeiDVmdot, VeiDVBnorm, VeiDVmax, VeiDVBasum, VeiDVBdot,
              VeiDVmdot, VeiDVBscal, VeiDVBcopy, VeiDVset, VeiDVBswap,
              VeiDVBaxpy, VeiDVmaxpy, VeiDVaypx, VeiDVwaxpy, VeiDVpmult,
              VeiDVpdiv,  
              VeiDVinsertvalues,0,0,
              VeiDVgetarray, VeiDVsize, VeiDVsize,VeiDVrange };

int VecCreateSequentialBLAS(int n,Vec *V)
{
  int      size = sizeof(DvVector)+n*sizeof(Scalar);
  Vec      v;
  DvVector *s;
  *V             = 0;
  PETSCHEADERCREATE(v,_Vec,VEC_COOKIE,SEQVECTOR,MPI_COMM_SELF);
  PLogObjectCreate(v);
  v->destroy     = VeiDestroyVector;
  v->view        = VeiDVview;
  s              = (DvVector *) MALLOC(size); CHKPTR(s);
  v->ops         = &DvOps;
  v->data        = (void *) s;
  s->n           = n;
  s->array       = (Scalar *)(s + 1);
  *V = v; return 0;
}

static int VeiDVBCreateVector(Vec win,Vec *V)
{
  DvVector *w = (DvVector *)win->data;
  return VecCreateSequentialBLAS(w->n,V);
}

