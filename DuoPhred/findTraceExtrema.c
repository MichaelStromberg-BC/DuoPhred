/** findTraceExtrema.c **/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "phred.h"

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

#ifdef ANSI_C
int findTraceExtrema( ChromatData *chromatData )
#else
int findTraceExtrema( chromatData )
ChromatData *chromatData;
#endif
{
  int i, j;
  int numPoint;
  FLOAT min, max;
  FLOAT *trace;

  min      = chromatData->trace[0][0];
  max      = chromatData->trace[0][0];
  numPoint = chromatData->numPoint;

  for( j = 0; j < 4; ++j )
  {
    trace = chromatData->trace[j];
    for( i = 0; i < numPoint; ++i )
    {
      if( trace[i] < min )
      {
        min = trace[i];
      }

      if( trace[i] > max )
      {
        max = trace[i];
      }
    }
  }

  chromatData->minTraceValue = min;
  chromatData->maxTraceValue = max;

  return( OK );
}
