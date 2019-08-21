static char help[]= "Test PetscSFFCompose when the ilocal array is not the identity\n\n";

#include <petsc.h>
#include <petscsf.h>

int main(int argc, char **argv)
{
  PetscErrorCode     ierr;
  PetscSF            sfA, sfB, sfBA;
  const PetscInt     nroots = 4, nleaves = 4;
  PetscInt          *ilocal;
  PetscSFNode       *iremote;
  Vec                a, b, ba;
  const PetscScalar *arrayR;
  PetscScalar       *arrayW;
  PetscMPIInt        rank, size;

  ierr = PetscInitialize(&argc,&argv,NULL,help);if (ierr) return ierr;
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
  ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);

  if (size != 1) SETERRQ(PETSC_COMM_WORLD, PETSC_ERR_USER, "Only coded for one MPI process");

  ierr = PetscSFCreate(PETSC_COMM_WORLD, &sfA); CHKERRQ(ierr);

  ierr = PetscSFCreate(PETSC_COMM_WORLD, &sfB); CHKERRQ(ierr);

  ierr = PetscMalloc1(nleaves, &ilocal); CHKERRQ(ierr);
  ierr = PetscMalloc1(nleaves, &iremote); CHKERRQ(ierr);

  for ( PetscInt i = 0; i < nleaves; i++ ) {
    ilocal[i] = nleaves - i - 1;
    iremote[i].rank = 0;
    iremote[i].index = i;
  }

  ierr = PetscSFSetGraph(sfA, nroots, nleaves, ilocal, PETSC_COPY_VALUES, iremote, PETSC_COPY_VALUES); CHKERRQ(ierr);

  for ( PetscInt i = 0; i < nleaves; i++ ) {
    ilocal[i] = i;
    iremote[i].rank = 0;
    iremote[i].index = i;
  }

  ierr = PetscSFSetGraph(sfB, nroots, nleaves, ilocal, PETSC_OWN_POINTER, iremote, PETSC_OWN_POINTER); CHKERRQ(ierr);

  ierr = PetscSFSetUp(sfA); CHKERRQ(ierr);

  ierr = PetscSFSetUp(sfB); CHKERRQ(ierr);

  ierr = PetscObjectSetName((PetscObject)sfA, "sfA");CHKERRQ(ierr);
  ierr = PetscObjectSetName((PetscObject)sfB, "sfB");CHKERRQ(ierr);

  ierr = VecCreateSeq(PETSC_COMM_WORLD, nroots, &a); CHKERRQ(ierr);
  ierr = VecCreateSeq(PETSC_COMM_WORLD, nleaves, &b); CHKERRQ(ierr);
  ierr = VecDuplicate(b, &ba); CHKERRQ(ierr);
  ierr = VecGetArray(a, &arrayW); CHKERRQ(ierr);
  for ( PetscInt i = 0; i < nroots; i++ ) {
    arrayW[i] = (PetscScalar)i;
  }
  ierr = VecRestoreArray(a, &arrayW); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "Initial Vec A\n"); CHKERRQ(ierr);
  ierr = VecView(a, NULL); CHKERRQ(ierr);
  ierr = VecGetArrayRead(a, &arrayR); CHKERRQ(ierr);
  ierr = VecGetArray(b, &arrayW); CHKERRQ(ierr);

  ierr = PetscSFBcastBegin(sfA, MPIU_SCALAR, arrayR, arrayW); CHKERRQ(ierr);
  ierr = PetscSFBcastEnd(sfA, MPIU_SCALAR, arrayR, arrayW); CHKERRQ(ierr);
  ierr = VecRestoreArray(b, &arrayW); CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(a, &arrayR); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "\nBroadcast A->B over sfA\n"); CHKERRQ(ierr);
  ierr = VecView(b, NULL); CHKERRQ(ierr);

  ierr = VecGetArrayRead(b, &arrayR); CHKERRQ(ierr);
  ierr = VecGetArray(ba, &arrayW); CHKERRQ(ierr);
  ierr = PetscSFBcastBegin(sfB, MPIU_SCALAR, arrayR, arrayW); CHKERRQ(ierr);
  ierr = PetscSFBcastEnd(sfB, MPIU_SCALAR, arrayR, arrayW); CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(b, &arrayR); CHKERRQ(ierr);
  ierr = VecRestoreArray(ba, &arrayW); CHKERRQ(ierr);

  ierr = PetscPrintf(PETSC_COMM_WORLD, "\nBroadcast B->BA over sfB\n"); CHKERRQ(ierr);
  ierr = VecView(ba, NULL); CHKERRQ(ierr);

  ierr = PetscSFCompose(sfA, sfB, &sfBA); CHKERRQ(ierr);

  ierr = PetscObjectSetName((PetscObject)sfBA, "(sfB o sfA)");CHKERRQ(ierr);
  ierr = VecGetArrayRead(a, &arrayR); CHKERRQ(ierr);
  ierr = VecGetArray(ba, &arrayW); CHKERRQ(ierr);

  ierr = PetscSFBcastBegin(sfBA, MPIU_SCALAR, arrayR, arrayW); CHKERRQ(ierr);
  ierr = PetscSFBcastEnd(sfBA, MPIU_SCALAR, arrayR, arrayW); CHKERRQ(ierr);
  ierr = VecRestoreArray(ba, &arrayW); CHKERRQ(ierr);
  ierr = VecRestoreArrayRead(a, &arrayR); CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD, "\nBroadcast A->BA over sfBA (sfB o sfA)\n"); CHKERRQ(ierr);
  ierr = VecView(ba, NULL); CHKERRQ(ierr);

  ierr = VecDestroy(&ba); CHKERRQ(ierr);
  ierr = VecDestroy(&b); CHKERRQ(ierr);
  ierr = VecDestroy(&a); CHKERRQ(ierr);

  ierr = PetscSFView(sfA, NULL); CHKERRQ(ierr);
  ierr = PetscSFView(sfB, NULL); CHKERRQ(ierr);
  ierr = PetscSFView(sfBA, NULL); CHKERRQ(ierr);
  ierr = PetscSFDestroy(&sfA); CHKERRQ(ierr);
  ierr = PetscSFDestroy(&sfB); CHKERRQ(ierr);
  ierr = PetscSFDestroy(&sfBA); CHKERRQ(ierr);

  ierr = PetscFinalize();

  return ierr;
}

/*TEST

   test:

TEST*/

