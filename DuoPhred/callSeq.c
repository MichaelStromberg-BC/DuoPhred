/** callSeq.c **/

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


/*
** Call bases.
*/
#ifdef ANSI_C
int callSeq( char *inFileName, PhredData *phredData, Option *option, int *status )
#else
int callSeq( inFileName, phredData, option, status )
char *inFileName;
PhredData *phredData;
Option *option;
int *status;
#endif
{
  int i;
  int istat;
  int numByte;
  int numPoint;
  register double scaleFactor;
  FLOAT **tr_vals;
  char line[2*PHRED_PATH_MAX];

  numPoint = phredData->numPoint;

  *status = 0;

  /*
  ** Are there at least a minimum number of
  ** trace points?
  */
  if( numPoint < 500 )
  {
    if( option->tagOption )
    {
/*
** OK
*/
      sprintf( line,
               "PROCESSING_ERROR: %.*s: unable to call bases: trace length %d < 500 point minimum\n",
               PHRED_PATH_MAX,
               inFileName,
               phredData->numPoint );
    }
    else
    {
/*
** OK
*/
      sprintf( line, "unable to call bases: trace length %d < 500 point minimum\n", phredData->numPoint );
    }
    fprintf( stderr, "%s", line );
    if( option->log )
    {
      writeLog( line );
    }
    phredData->leftTrimPoint = 0;
    phredData->rghtTrimPoint = 0;
    *status = 1;
    return( ERROR );
  }

  /*
  ** Allocate memory for scaled trace data.
  */
  numByte = 4 * sizeof( FLOAT * );
  tr_vals     = (FLOAT **)ourMalloc( numByte );

  numByte = numPoint * sizeof( FLOAT );
  for( i = 0; i < 4; ++i )
  {
    tr_vals[i] = (FLOAT *)ourMalloc( numByte );
  }

  scaleFactor = MAX_TRACE_VAL / phredData->maxTraceValue;

  for( i = 0; i < numPoint; ++i )
  {
    tr_vals[0][i] = phredData->trace[0][i] * scaleFactor;
    tr_vals[1][i] = phredData->trace[1][i] * scaleFactor;
    tr_vals[2][i] = phredData->trace[2][i] * scaleFactor;
    tr_vals[3][i] = phredData->trace[3][i] * scaleFactor;
  }


  fit_main( numPoint, tr_vals, phredData, &istat );
  if( istat == ERROR )
  {
    for( i = 0; i < 4; ++i )
    {
      ourFree( (char *)tr_vals[i] );
    }
    ourFree( (char *)tr_vals );
    phredData->leftTrimPoint = 0;
    phredData->rghtTrimPoint = 0;
    *status = 1;
    return( ERROR );
  }

  /*
  ** free up the dynamic memory
  */
  for( i = 0; i < 4; ++i )
  {
    ourFree( (char *)tr_vals[i] );
  }
  ourFree( (char *)tr_vals );

  *status = 0;

  return( OK );

}

