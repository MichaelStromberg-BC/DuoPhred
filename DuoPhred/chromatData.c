/** chromatData.c **/

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
#include "ourMalloc.h"
#include "ourFree.h"

/*
** ChromatData allocate/free.
*/
#ifdef ANSI_C
ChromatData *allocChromatData( int numPoint, int numBase )
#else
ChromatData *allocChromatData( numPoint, numBase )
int numPoint;
int numBase;
#endif
{
  int j;
  int numByte;
  ChromatData *chromatData;

  numByte = sizeof( ChromatData );
  chromatData = (ChromatData *)ourMalloc( numByte );

  for( j = 0; j < 4; ++j )
  {
    chromatData->trace[j] = NULL;
  }
  chromatData->base     = NULL;
  chromatData->baseLoc  = NULL;
  chromatData->baseQual = NULL;

  if( numPoint > 0 )
  {
    numByte = numPoint * sizeof( FLOAT );
  }
  else
  {
    numByte = sizeof( FLOAT );
  }

  for( j = 0; j < 4; ++j )
  {
    chromatData->trace[j] = (FLOAT *)ourMalloc( numByte );
  }

  if( numBase > 0 )
  {
    numByte = numBase * sizeof( char );
  }
  else
  {
    numByte = sizeof( char );
  }
  chromatData->base = (char *)ourMalloc( numByte );

  if( numBase > 0 )
  {
    numByte = numBase * sizeof( int );
  }
  else
  {
    numByte = sizeof( int );
  }
  chromatData->baseLoc = (int *)ourMalloc( numByte );

  if( numBase > 0 )
  {
    numByte = numBase * sizeof( int );
  }
  else
  {
    numByte = sizeof( int );
  }
  chromatData->baseQual = (int *)ourMalloc( numByte );

  /*
  ** Initialize values.
  */
  chromatData->fileName[0]    = '\0';
  chromatData->fileType       = -1;
  chromatData->numPoint       = numPoint;
  chromatData->maxTraceValue  = 0.0;
  chromatData->numBase        = numBase;
  chromatData->primerLoc      = 0;
  chromatData->avgSpacing     = 0.0;
  chromatData->sampleName[0]  = '\0';
  chromatData->machineName[0] = '\0';
  chromatData->primerID[0]    = '\0'; 

  for( j = 0; j < 4; ++j )
  {
    chromatData->signalStrength[j] = 0;
  }

  return( chromatData );
}

#ifdef ANSI_C
int freeChromatData( ChromatData *chromatData )
#else
int freeChromatData( chromatData )
ChromatData *chromatData;
#endif
{
  int j;

  for( j = 0; j < 4; ++j )
  {
    if( chromatData->trace[j] != NULL )
    {
      ourFree( (char *)chromatData->trace[j] );
    }
  }

  if( chromatData->base != NULL )
  {
    ourFree( (char *)chromatData->base );
  }

  if( chromatData->baseLoc != NULL )
  {
    ourFree( (char *)chromatData->baseLoc );
  }

  if( chromatData->baseQual != NULL )
  {
    ourFree( (char *)chromatData->baseQual );
  }

  ourFree( (char *)chromatData );

  return( OK );
}

