#ifdef PETSC_RCS_HEADER
static char vcid[] = "$Id: inherit.c,v 1.52 1999/07/15 14:13:35 bsmith Exp bsmith $";
#endif
/*
     Provides utility routines for manipulating any type of PETSc object.
*/
#include "petsc.h"  /*I   "petsc.h"    I*/

extern int PetscObjectGetComm_Petsc(PetscObject,MPI_Comm *);
extern int PetscObjectCompose_Petsc(PetscObject,const char[],PetscObject);
extern int PetscObjectQuery_Petsc(PetscObject,const char[],PetscObject *);
extern int PetscObjectComposeFunction_Petsc(PetscObject,const char[],const char[],void *);
extern int PetscObjectQueryFunction_Petsc(PetscObject,const char[],void **);
extern int PetscObjectComposeLanguage_Petsc(PetscObject,PetscLanguage,void *);
extern int PetscObjectQueryLanguage_Petsc(PetscObject,PetscLanguage,void **);

#undef __FUNC__  
#define __FUNC__ "PetscHeaderCreate_Private"
/*
   PetscHeaderCreate_Private - Creates a base PETSc object header and fills
   in the default values.  Called by the macro PetscHeaderCreate().
*/
int PetscHeaderCreate_Private(PetscObject h,int cookie,int type,char *class_name,MPI_Comm comm,
                              int (*des)(PetscObject),int (*vie)(PetscObject,Viewer))
{
  static int idcnt = 1;
  int        ierr;

  PetscFunctionBegin;
  h->cookie                 = cookie;
  h->type                   = type;
  h->class_name             = class_name;
  h->prefix                 = 0;
  h->refct                  = 1;
  h->amem                   = -1;
  h->id                     = idcnt++;
  h->parentid               = 0;
  h->bops->destroy          = des;
  h->bops->view             = vie;
  h->bops->getcomm          = PetscObjectGetComm_Petsc;
  h->bops->compose          = PetscObjectCompose_Petsc;
  h->bops->query            = PetscObjectQuery_Petsc;
  h->bops->composefunction  = PetscObjectComposeFunction_Petsc;
  h->bops->queryfunction    = PetscObjectQueryFunction_Petsc;
  h->bops->querylanguage    = PetscObjectQueryLanguage_Petsc;
  h->bops->composelanguage  = PetscObjectComposeLanguage_Petsc;
  ierr = PetscCommDuplicate_Private(comm,&h->comm,&h->tag);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "PetscHeaderDestroy_Private"
/*
    PetscHeaderDestroy_Private - Destroys a base PETSc object header. Called by 
    the macro PetscHeaderDestroy().
*/
int PetscHeaderDestroy_Private(PetscObject h)
{
  int ierr;

  PetscFunctionBegin;
  if (h->amem != -1) {
    SETERRQ(1,1,"PETSc object destroyed before its AMS publication was destroyed");
  }

  ierr = PetscCommDestroy_Private(&h->comm);CHKERRQ(ierr);
  ierr = PetscFree(h->bops);CHKERRQ(ierr);
  ierr = PetscFree(h->ops);CHKERRQ(ierr);
  ierr = OListDestroy(&h->olist);CHKERRQ(ierr);
  ierr = FListDestroy(h->qlist);CHKERRQ(ierr);
  if (h->type_name) {ierr = PetscFree(h->type_name);CHKERRQ(ierr);}
  if (h->name) {ierr = PetscFree(h->name);CHKERRQ(ierr);}
  h->cookie = PETSCFREEDHEADER;
  if (h->prefix) {ierr = PetscFree(h->prefix);CHKERRQ(ierr);}
  if (h->fortran_func_pointers) {
    ierr = PetscFree(h->fortran_func_pointers);CHKERRQ(ierr);
  }
  ierr = PetscFree(h);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "PetscObjectReference"
/*@C
   PetscObjectReference - Indicates to any PetscObject that it is being
   referenced by another PetscObject. This increases the reference
   count for that object by one.

   Collective on PetscObject

   Input Parameter:
.  obj - the PETSc object. This must be cast with (PetscObject), for example, 
         PetscObjectReference((PetscObject) mat);

   Level: advanced

.seealso: PetscObjectCompose(), PetscObjectDereference()
@*/
int PetscObjectReference(PetscObject obj)
{
  PetscFunctionBegin;
  PetscValidHeader(obj);
  obj->refct++;
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "PetscObjectGetReference"
/*@C
   PetscObjectGetReference - Gets the current reference count for 
   any PETSc object.

   Not Collective

   Input Parameter:
.  obj - the PETSc object; this must be cast with (PetscObject), for example, 
         PetscObjectGetReference((PetscObject) mat,&cnt);

   Output Parameter:
.  cnt - the reference count

   Level: advanced

.seealso: PetscObjectCompose(), PetscObjectDereference(), PetscObjectReference()
@*/
int PetscObjectGetReference(PetscObject obj,int *cnt)
{
  PetscFunctionBegin;
  PetscValidHeader(obj);
  *cnt = obj->refct;
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "PetscObjectDereference"
/*@
   PetscObjectDereference - Indicates to any PetscObject that it is being
   referenced by one less PetscObject. This decreases the reference
   count for that object by one.

   Collective on PetscObject

   Input Parameter:
.  obj - the PETSc object; this must be cast with (PetscObject), for example, 
         PetscObjectDereference((PetscObject) mat);

   Level: advanced

.seealso: PetscObjectCompose(), PetscObjectReference()
@*/
int PetscObjectDereference(PetscObject obj)
{
  int ierr;

  PetscFunctionBegin;
  PetscValidHeader(obj);
  if (obj->bops->destroy) {
    ierr = (*obj->bops->destroy)(obj);CHKERRQ(ierr);
  } else if (--obj->refct == 0) {
    SETERRQ(PETSC_ERR_SUP,0,"This PETSc object does not have a generic destroy routine");
  }
  PetscFunctionReturn(0);
}

/* ----------------------------------------------------------------------- */
/*
     The following routines are the versions private to the PETSc object
     data structures.
*/
#undef __FUNC__  
#define __FUNC__ "PetscObjectGetComm_Petsc"
int PetscObjectGetComm_Petsc(PetscObject obj,MPI_Comm *comm)
{
  PetscFunctionBegin;
  *comm = obj->comm;
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "PetscObjectCompose_Petsc"
int PetscObjectCompose_Petsc(PetscObject obj,const char name[],PetscObject ptr)
{
  int  ierr;
  char *tname;

  PetscFunctionBegin;
  if (ptr) {
    ierr = OListReverseFind(ptr->olist,obj,&tname);CHKERRQ(ierr);
    if (tname){
      SETERRQ(1,1,"An object cannot be composed with an object that was compose with it");
    }
  }
  ierr = OListAdd(&obj->olist,name,ptr);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "PetscObjectQuery_Petsc"
int PetscObjectQuery_Petsc(PetscObject obj,const char name[],PetscObject *ptr)
{
  int ierr;

  PetscFunctionBegin;
  ierr = OListFind(obj->olist,name,ptr);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "PetscObjectComposeLanguage_Petsc"
int PetscObjectComposeLanguage_Petsc(PetscObject obj,PetscLanguage lang,void *vob)
{
  PetscFunctionBegin;
  if (lang == PETSC_LANGUAGE_CPP) {
    obj->cpp = vob;
  } else {
    SETERRQ(1,1,"No support for this language yet");
  }
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "PetscObjectQueryLanguage_Petsc"
int PetscObjectQueryLanguage_Petsc(PetscObject obj,PetscLanguage lang,void **vob)
{
  PetscFunctionBegin;
  if (lang == PETSC_LANGUAGE_C) {
    *vob = (void *) obj;
  } else if (lang == PETSC_LANGUAGE_CPP) {
    if (obj->cpp) {
      *vob = obj->cpp;
    } else {
      SETERRQ(1,1,"No C++ wrapper generated");
    }
  } else {
    SETERRQ(1,1,"No support for this language yet");
  }
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "PetscObjectComposeFunction_Petsc"
int PetscObjectComposeFunction_Petsc(PetscObject obj,const char name[],
                                     const char fname[],void *ptr)
{
  int ierr;

  PetscFunctionBegin;
  ierr = FListAdd(&obj->qlist,name,fname,(int (*)(void *))ptr);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "PetscObjectQueryFunction_Petsc"
int PetscObjectQueryFunction_Petsc(PetscObject obj,const char name[],void **ptr)
{
  int ierr;

  PetscFunctionBegin;
  ierr = FListFind(obj->comm,obj->qlist,name,( int(**)(void *)) ptr);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/*
        These are the versions that are usable to any CCA compliant objects
*/
#undef __FUNC__  
#define __FUNC__ "PetscObjectCompose"
/*@C
   PetscObjectCompose - Associates another PETSc object with a given PETSc object. 
                       
   Not Collective

   Input Parameters:
+  obj - the PETSc object; this must be cast with (PetscObject), for example, 
         PetscObjectCompose((PetscObject) mat,...);
.  name - name associated with the child object 
-  ptr - the other PETSc object to associate with the PETSc object; this must also be 
         cast with (PetscObject)

   Level: advanced

   Notes:
   The second object's reference count is automatically increased by one when it is
   composed.

   Replaces any previous object that had the same name.

   If ptr is null and name has previously been composed using an object, then that
   entry is removed from the obj.

   PetscObjectCompose() can be used with any PETSc object (such as
   Mat, Vec, KSP, SNES, etc.) or any user-provided object.  See 
   PetscObjectContainerCreate() for info on how to create an object from a 
   user-provided pointer that may then be composed with PETSc objects.
   
.keywords: object, composition

.seealso: PetscObjectQuery(), PetscObjectContainerCreate()
@*/
int PetscObjectCompose(PetscObject obj,const char name[],PetscObject ptr)
{
  int ierr;

  PetscFunctionBegin;
  ierr = (*obj->bops->compose)(obj,name,ptr);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "PetscObjectQuery"
/*@C
   PetscObjectQuery  - Gets a PETSc object associated with a given object.
                       
   Not Collective

   Input Parameters:
+  obj - the PETSc object
         Thus must be cast with a (PetscObject), for example, 
         PetscObjectCompose((PetscObject) mat,...);
.  name - name associated with child object 
-  ptr - the other PETSc object associated with the PETSc object, this must also be 
         cast with (PetscObject)

   Level: advanced

.keywords: object, composition

.seealso: PetscObjectQuery()
@*/
int PetscObjectQuery(PetscObject obj,const char name[],PetscObject *ptr)
{
  int ierr;

  PetscFunctionBegin;
  ierr = (*obj->bops->query)(obj,name,ptr);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "PetscObjectQueryLanguage"
/*@C
   PetscObjectQueryLanguage - Returns a language specific interface to the given object
                       
   Not Collective

   Input Parameters:
+  obj - the PETSc object
         Thus must be cast with a (PetscObject), for example, 
         PetscObjectCompose((PetscObject) mat,...);
-  lang - one of PETSC_LANGUAGE_C, PETSC_LANGUAGE_F77, PETSC_LANGUAGE_CPP

   Output Parameter:
.  ptr - the language specific interface

   Level: developer

.keywords: language, object wrappers

.seealso: PetscObjectQuery()
@*/
int PetscObjectQueryLanguage(PetscObject obj,PetscLanguage lang,void **ptr)
{
  int ierr;

  PetscFunctionBegin;
  ierr = (*obj->bops->querylanguage)(obj,lang,ptr);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "PetscObjectComposeLanguage"
/*@C
   PetscObjectComposeLanguage - Sets a language specific interface to the given object
                       
   Not Collective

   Input Parameters:
+  obj - the PETSc object
         Thus must be cast with a (PetscObject), for example, 
         PetscObjectCompose((PetscObject) mat,...);
.  lang - one of PETSC_LANGUAGE_C, PETSC_LANGUAGE_F77, PETSC_LANGUAGE_CPP
-  ptr - the language specific interface

   Level: developer

.keywords: language, object wrappers

.seealso: PetscObjectQuery()
@*/
int PetscObjectComposeLanguage(PetscObject obj,PetscLanguage lang,void *ptr)
{
  int ierr;

  PetscFunctionBegin;
  ierr = (*obj->bops->composelanguage)(obj,lang,ptr);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}


#undef __FUNC__  
#define __FUNC__ "PetscObjectComposeFunction_Private"
int PetscObjectComposeFunction_Private(PetscObject obj,const char name[],
                                       const char fname[],void *ptr)
{
  int ierr;

  PetscFunctionBegin;
  ierr = (*obj->bops->composefunction)(obj,name,fname,ptr);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "PetscObjectQueryFunction"
/*@C
   PetscObjectQueryFunction - Gets a function associated with a given object.
                       
   Collective on PetscObject

   Input Parameters:
+  obj - the PETSc object; this must be cast with (PetscObject), for example, 
         PetscObjectQueryFunction((PetscObject) ksp,...);
-  name - name associated with the child function

   Output Parameter:
.  ptr - function pointer

   Level: advanced

.keywords: object, composition

.seealso: PetscObjectComposeFunction()
@*/
int PetscObjectQueryFunction(PetscObject obj,const char name[],void **ptr)
{
  int ierr;

  PetscFunctionBegin;
  ierr = (*obj->bops->queryfunction)(obj,name,ptr);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

/* -------------------------------------------------------------------------------------*/

#undef __FUNC__  
#define __FUNC__ "PetscDataTypeToMPIDataType"
int PetscDataTypeToMPIDataType(PetscDataType ptype,MPI_Datatype* mtype)
{
  PetscFunctionBegin;
  if (ptype == PETSC_INT) {
    *mtype = MPI_INT;
  } else if (ptype == PETSC_DOUBLE) {
    *mtype = MPI_DOUBLE;
  } else if (ptype == PETSC_SCALAR) {
    *mtype = MPIU_SCALAR;
#if defined(PETSC_USE_COMPLEX)
  } else if (ptype == PETSC_COMPLEX) {
    *mtype = MPIU_COMPLEX;
#endif
  } else if (ptype == PETSC_CHAR) {
    *mtype = MPI_CHAR;
  } else if (ptype == PETSC_LOGICAL) {
    *mtype = MPI_BYTE;
  } else {
    SETERRQ(PETSC_ERR_ARG_OUTOFRANGE,1,"Unknown PETSc datatype");
  }
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "PetscDataTypeGetSize"
int PetscDataTypeGetSize(PetscDataType ptype,int *size)
{
  PetscFunctionBegin;
  if (ptype == PETSC_INT) {
    *size = PETSC_INT_SIZE;
  } else if (ptype == PETSC_DOUBLE) {
    *size = PETSC_DOUBLE_SIZE;
  } else if (ptype == PETSC_SCALAR) {
    *size = PETSC_SCALAR_SIZE;
#if defined(PETSC_USE_COMPLEX)
  } else if (ptype == PETSC_COMPLEX) {
    *size = PETSC_COMPLEX_SIZE;
#endif
  } else if (ptype == PETSC_CHAR) {
    *size = PETSC_CHAR_SIZE;
  } else if (ptype == PETSC_LOGICAL) {
    *size = PETSC_LOGICAL_SIZE;
  } else {
    SETERRQ(PETSC_ERR_ARG_OUTOFRANGE,1,"Unknown PETSc datatype");
  }
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "PetscDataTypeGetName"
int PetscDataTypeGetName(PetscDataType ptype,char *name[])
{
  PetscFunctionBegin;
  if (ptype == PETSC_INT) {
    *name = "int";
  } else if (ptype == PETSC_DOUBLE) {
    *name = "double";
  } else if (ptype == PETSC_SCALAR) {
    *name = "Scalar";
#if defined(PETSC_USE_COMPLEX)
  } else if (ptype == PETSC_COMPLEX) {
    *name = "complex";
#endif
  } else if (ptype == PETSC_CHAR) {
    *name = "char";
  } else if (ptype == PETSC_LOGICAL) {
    *name = "logical";
  } else {
    SETERRQ(PETSC_ERR_ARG_OUTOFRANGE,1,"Unknown PETSc datatype");
  }
  PetscFunctionReturn(0);
}

struct _p_PetscObjectContainer {
  PETSCHEADER(int)
  void   *ptr;
};

#undef __FUNC__  
#define __FUNC__ "PetscObjectContainerGetPointer"
/*@C
   PetscObjectContainerGetPointer - Gets the pointer value contained in the container.

   Collective on PetscObjectContainer

   Input Parameter:
.  obj - the object created with PetscObjectContainerCreate()

   Output Parameter:
.  ptr - the pointer value

   Level: advanced

.seealso: PetscObjectContainerCreate(), PetscObjectContainerDestroy(), 
          PetscObjectContainerSetPointer()
@*/
int PetscObjectContainerGetPointer(PetscObjectContainer obj,void **ptr)
{
  PetscFunctionBegin;
  *ptr = obj->ptr;
  PetscFunctionReturn(0);
}


#undef __FUNC__  
#define __FUNC__ "PetscObjectContainerSetPointer"
/*@C
   PetscObjectContainerSetPointer - Sets the pointer value contained in the container.

   Collective on PetscObjectContainer

   Input Parameters:
+  obj - the object created with PetscObjectContainerCreate()
-  ptr - the pointer value

   Level: advanced

.seealso: PetscObjectContainerCreate(), PetscObjectContainerDestroy(), 
          PetscObjectContainerGetPointer()
@*/
int PetscObjectContainerSetPointer(PetscObjectContainer obj,void *ptr)
{
  PetscFunctionBegin;
  obj->ptr = ptr;
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "PetscObjectContainerDestroy"
/*@C
   PetscObjectContainerDestroy - Destroys a PETSc container object.

   Collective on PetscObjectContainer

   Input Parameter:
.  obj - an object that was created with PetscObjectContainerCreate()

   Level: advanced

.seealso: PetscObjectContainerCreate()
@*/
int PetscObjectContainerDestroy(PetscObjectContainer obj)
{
  PetscFunctionBegin;
  if (--obj->refct > 0) PetscFunctionReturn(0);
  PetscHeaderDestroy(obj);
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ "PetscObjectContainerCreate"
/*@C
   PetscObjectContainerCreate - Creates a PETSc object that has room to hold
   a single pointer. This allows one to attach any type of data (accessible
   through a pointer) with the PetscObjectCompose() function to a PetscObject.

   Collective on MPI_Comm

   Input Parameters:
.  comm - MPI communicator that shares the object

   Level: advanced

.seealso: PetscObjectContainerDestroy()
@*/
int PetscObjectContainerCreate(MPI_Comm comm,PetscObjectContainer *container)
{
  PetscObjectContainer contain;

  PetscFunctionBegin;

  PetscHeaderCreate(contain, _p_PetscObjectContainer,int,PETSC_COOKIE,0,"container",
                    comm,PetscObjectContainerDestroy,0);
  *container = contain;
  PetscFunctionReturn(0);
}

/*MC
   PetscObjectComposeFunction - Associates a function with a given PETSc object. 
                       
   Collective on PetscObject

   Input Parameters:
+  obj - the PETSc object; this must be cast with a (PetscObject), for example, 
         PetscObjectCompose((PetscObject) mat,...);
.  name - name associated with the child function
.  fname - name of the function
-  ptr - function pointer (or PETSC_NULL if using dynamic libraries)

   Level: advanced

    Synopsis:
    PetscObjectComposeFunction(PetscObject obj,char *name,char *fname,void *ptr)

   Notes:
   PetscObjectComposeFunction() can be used with any PETSc object (such as
   Mat, Vec, KSP, SNES, etc.) or any user-provided object. 

   The composed function must be wrapped in a EXTERN_C_BEGIN/END for this to
   work in C++/complex with dynamic link libraries (PETSC_USE_DYNAMIC_LIBRARIES)
   enabled.

.keywords: object, composition

.seealso: PetscObjectQueryFunction()
M*/
