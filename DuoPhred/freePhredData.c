/** freePhredData.c **/

/*
*|***************************************************************************|*
*|                                                                           |*
*|   Program: phred                                                          |*
*|   Version: 0.000925.c                                                     |*
*|                                                                           |*
*|   Copyright (C) 1993-2000 by Phil Green and Brent Ewing.                  |*
*|   All rights reserved.                                                    |*
*|                                                                           |*
*|   This software is a beta-test version of the phred package.              |*
*|   It should not be redistributed or used for any commercial               |*
*|   purpose, including commercially funded sequencing, without              |*
*|   written permission from the author and the University of                |*
*|   Washington.                                                             |*
*|                                                                           |*
*|   This software is provided ``AS IS'' and any express or                  |*
*|   implied warranties, including, but not limited to, the                  |*
*|   implied warranties of merchantability and fitness for a                 |*
*|   particular purpose, are disclaimed.  In no event shall                  |*
*|   the authors or the University of Washington be liable for               |*
*|   any direct, indirect, incidental, special, exemplary, or                |*
*|   consequential damages (including, but not limited to,                   |*
*|   procurement of substitute goods or services; loss of use,               |*
*|   data, or profits; or business interruption) however caused              |*
*|   and on any theory of liability, whether in contract, strict             |*
*|   liability, or tort (including negligence or otherwise)                  |*
*|   arising in any way out of the use of this software, even                |*
*|   if advised of the possibility of such damage.                           |*
*|                                                                           |*
*|   Portions of the code benefit from ideas due to Dave Ficenec,            |*
*|   LaDeana Hillier, Mike Wendl, and Tim Gleeson.  These are                |*
*|   indicated in the relevant source files.                                 |*
*|                                                                           |*
*|***************************************************************************|*
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "phred.h"

#ifdef ANSI_C
int freePhredData( PhredData *phredData )
#else
int freePhredData( phredData )
PhredData *phredData;
#endif
{
  int j;

  for( j = 0; j < 4; ++j )
  {
    if( phredData->trace[j] != NULL )
    {
      ourFree( (char *)phredData->trace[j] );
    }
  }

  for( j = 0; j < 2; ++j )
  {
    if( phredData->base[j] != NULL )
    {
      ourFree( (char *)phredData->base[j] );
    }
    if( phredData->baseLoc[j] != NULL )
    {
      ourFree( (char *)phredData->baseLoc[j] );
    }
    if( phredData->baseQual[j] != NULL )
    {
      ourFree( (char *)phredData->baseQual[j] );
    }
  }

  if( phredData->baseDsc != NULL )
  {
    ourFree( (char *)phredData->baseDsc );
  }
  if( phredData->qualIndex != NULL )
  {
    ourFree( (char *)phredData->qualIndex );
  }

  ourFree( (char *)phredData );

  return( OK );
}
