#ifndef lint
static char vcid[] = "$Id: dsetpause.c,v 1.6 1996/08/08 14:44:45 bsmith Exp balay $";
#endif
/*
       Provides the calling sequences for all the basic Draw routines.
*/
#include "src/draw/drawimpl.h"  /*I "draw.h" I*/

#undef __FUNCTION__  
#define __FUNCTION__ "DrawSetPause"
/*@
   DrawSetPause - Sets the amount of time that program pauses after 
   a DrawPause() is called. 

   Input Paramters:
.  draw - the drawing object
.  pause - number of seconds to pause, -1 implies until user input

   Note:
   By default the pause time is zero unless the -draw_pause option is given 
   during DrawOpenX().

.keywords: draw, set, pause

.seealso: DrawGetPause(), DrawPause()
@*/
int DrawSetPause(Draw draw,int pause)
{
  PetscValidHeaderSpecific(draw,DRAW_COOKIE);
  if (draw->type == DRAW_NULLWINDOW) return 0;
  draw->pause = pause;
  return 0;
}
