#ifdef PETSC_RCS_HEADER
static char vcid[] = "$Id: PetscTime.c,v 1.11 1999/03/19 21:24:35 bsmith Exp bsmith $";
#endif

#include "petsc.h"
#include "pinclude/ptime.h"

#undef __FUNC__
#define __FUNC__ "main"
int main( int argc, char **argv)
{
  PLogDouble x, y;
  int        i;

  PetscInitialize(&argc, &argv,0,0);
  /* To take care of paging effects */
  PetscTime(y);

  for ( i=0; i<2; i++ ) { 
    PetscTime(x);
    PetscTime(y); 
    PetscTime(y);
    PetscTime(y);
    PetscTime(y);
    PetscTime(y);
    PetscTime(y); 
    PetscTime(y);
    PetscTime(y);
    PetscTime(y);
    PetscTime(y);

    fprintf(stderr,"%-15s : %e sec\n","PetscTime",(y-x)/10.0);
  }
  PetscTime(x);
  ierr = PetscSleep(10);CHKERRA(ierr);
  PetscTime(y); 
  fprintf(stderr,"%-15s : %e sec - Slept for 10 sec \n","PetscTime",(y-x));

  PetscFinalize();
  PetscFunctionReturn(0);
}
