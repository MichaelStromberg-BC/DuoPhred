/** trimSet.c **/

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

/*
*******************************************************************************
**                                                                           **
**    * trimSet.c was written by LaDeana Hillier.                       *    **
**                                                                           **
*******************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "phred.h"

#define ABS(A)	( ( (A) < 0 ) ? -(A) : (A) )

#ifdef ANSI_C
int trimSet( PhredData *phredData )
#else
int trimSet( phredData )
PhredData *phredData;
#endif
{
  int i, k;
  int diff;
  int minDiff;
  int leftCutoffBasePos;
  int rightCutoffBasePos;

  if( phredData->rghtTrimPoint <= 0 )
  {
    phredData->rghtTrimPoint = 1;
  }

  leftCutoffBasePos  = phredData->baseLoc[FIL][phredData->leftTrimPoint];
  rightCutoffBasePos = phredData->baseLoc[FIL][phredData->numBase[FIL] -
                       phredData->rghtTrimPoint];


  /*
  ** Figure out what position at which the right
  ** cutoff has been set and then ask Phil's
  ** local base calling what number that corresponds
  ** to for him.
  */
  k = 0;
  minDiff = ABS( phredData->baseLoc[LCL][0] - rightCutoffBasePos );
  for( i = 0; i < phredData->numBase[LCL]; ++i )
  {
    diff = ABS( phredData->baseLoc[LCL][i] - rightCutoffBasePos );
    if( diff < minDiff )
    {
      minDiff = diff;
      k = i;
    }
  }
  phredData->rghtTrimPoint = phredData->numBase[LCL] - k;

  /*
  ** Figure out what position at which the left
  ** cutoff has been set and then ask Phil's
  ** local base calling what number that corresponds
  ** to for him.
  */
  k = 0;
  minDiff = ABS( phredData->baseLoc[LCL][0] - leftCutoffBasePos );
  for( i = 0; i < phredData->numBase[LCL]; ++i )
  {
    diff = ABS( phredData->baseLoc[LCL][i] - leftCutoffBasePos );
    if( diff < minDiff )
    {
      minDiff = diff;
      k = i;
    }
  }
  phredData->leftTrimPoint = k;

  /*
  ** Watch for pathological conditions.
  */
  if( phredData->baseLoc[LCL][phredData->numBase[LCL]-1] <= rightCutoffBasePos )
  {
    phredData->rghtTrimPoint = 0;
  }
  if( phredData->baseLoc[LCL][phredData->numBase[LCL]-1] <= leftCutoffBasePos )
  {
    phredData->leftTrimPoint = phredData->numBase[LCL] - 1;
  }

  /*
  ** This is a special case when the rightCutoffBasePos
  ** does not fall within the expected baseLoc[] locations.
  */
  if( phredData->rghtTrimPoint > phredData->numBase[LCL] )
  {
    phredData->rghtTrimPoint = phredData->numBase[LCL];
  }

  return( OK );
}

