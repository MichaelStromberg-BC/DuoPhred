/** readSCF3.c **/

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
**    * readSCF3.c                                                      *     **
**    * benefits from ideas in code written by LaDeana Hillier and     *     **
**    * Tim Gleeson.                                                   *     **
**                                                                           **
*******************************************************************************
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <sys/types.h>
#include "phred.h"
#include "rwUtil.h"
#include "chromatData.h"
#include "freeChromatData.h"
#include "readSCF.h"


#ifdef ANSI_C
static int readComment3( int numEntry, CommentEntry *ce, ChromatData *chromatData )
#else
static int readComment3( numEntry, ce, chromatData )
int numEntry;
CommentEntry *ce;
ChromatData *chromatData;
#endif
{
  int i;
  char *cptr;

  /*
  ** Copy desired comments.
  */
  for( i = 0; i < numEntry; ++i )
  {
    if( strcmp( ce[i].label, "PRIM" ) == 0 ||
        strcmp( ce[i].label, "primer_position" ) == 0 )
    {
      chromatData->primerLoc = atoi( ce[i].value );
      ce[i].used = 1;
      continue;
    }

    if( strcmp( ce[i].label, "SPAC" ) == 0 ||
        strcmp( ce[i].label, "avg_spacing" ) == 0 )
    {
      chromatData->avgSpacing = atof( ce[i].value );
      ce[i].used = 1;
      continue;
    }

    if( strcmp( ce[i].label, "MACH" ) == 0 ||
        strcmp( ce[i].label, "machine_name" ) == 0 )
    {
      pstrcpy( chromatData->machineName, ce[i].value, PHRED_MAX_STRING_LEN );
      ce[i].used = 1;
      continue;
    }

    if( strcmp( ce[i].label, "NAME" ) == 0 ||
        strcmp( ce[i].label, "sample_name" ) == 0 )
    {
      pstrcpy( chromatData->sampleName, ce[i].value, PHRED_MAX_STRING_LEN );
      ce[i].used = 1;
      continue;
    }

    if( strcmp( ce[i].label, "DYEP" ) == 0 ||
        strcmp( ce[i].label, "dye_primer" ) == 0 )
    {
      pstrcpy( chromatData->primerID, ce[i].value, PHRED_MAX_STRING_LEN );
      ce[i].used = 1;
      continue;
    }

    if( strcmp( ce[i].label, "SIGN" ) == 0 ||
        strcmp( ce[i].label, "avg_signal_strength" ) == 0 )
    {
      if( ( cptr = findSignal( ce[i].value, 'A' ) ) != NULL )
      {
        chromatData->signalStrength[0] = atoi( cptr );
      }

      if( ( cptr = findSignal( ce[i].value, 'C' ) ) != NULL )
      {
        chromatData->signalStrength[1] = atoi( cptr );
      }

      if( ( cptr = findSignal( ce[i].value, 'G' ) ) != NULL )
      {
        chromatData->signalStrength[2] = atoi( cptr );
      }

      if( ( cptr = findSignal( ce[i].value, 'T' ) ) != NULL )
      {
        chromatData->signalStrength[3] = atoi( cptr );
      }

      ce[i].used = 1;
      continue;
    }

    if( strcmp( ce[i].label, "GELN" ) == 0 ||
        strcmp( ce[i].label, "gel_name" ) == 0 )
    {
      pstrcpy( chromatData->gelName, ce[i].value, PHRED_MAX_STRING_LEN );
      ce[i].used = 1;
      continue;
    }

    if( strcmp( ce[i].label, "LANE" ) == 0 ||
        strcmp( ce[i].label, "lane_number" ) == 0 )
    {
      chromatData->laneNumber = atoi( ce[i].value );
      ce[i].used = 1;
      continue;
    }

    if( strcmp( ce[i].label, "PROC" ) == 0 ||
        strcmp( ce[i].label, "processing" ) == 0 )
    {
      pstrcpy( chromatData->processing, ce[i].value, PHRED_MAX_STRING_LEN );
      ce[i].used = 1;
      continue;
    }

    if( strcmp( ce[i].label, "RTRK" ) == 0 ||
        strcmp( ce[i].label, "retraker" ) == 0 )
    {
      pstrcpy( chromatData->reTracker, ce[i].value, PHRED_MAX_STRING_LEN );
      ce[i].used = 1;
      continue;
    }

    if( strcmp( ce[i].label, "COMM" ) == 0 ||
        strcmp( ce[i].label, "comment" ) == 0 )
    {
      pstrcpy( chromatData->comment, ce[i].value, PHRED_MAX_STRING_LEN );
      ce[i].used = 1;
      continue;
    }

    if( strcmp( ce[i].label, "CONV" ) == 0 ||
        strcmp( ce[i].label, "conversion_program" ) == 0 )
       
    {
      pstrcpy( chromatData->convProg, ce[i].value, PHRED_MAX_STRING_LEN );
      ce[i].used = 1;
      continue;
    }

    if( strcmp( ce[i].label, "SRCE" ) == 0 ||
        strcmp( ce[i].label, "source" ) == 0 )
    {
      pstrcpy( chromatData->source, ce[i].value, PHRED_MAX_STRING_LEN );
      ce[i].used = 1;
      continue;
    }
  }

  /*
  ** Append other comments to "comment".
  */
  if( chromatData->comment[strlen( chromatData->comment )-1] != '\n' )
  {
    if( strlen( chromatData->comment ) < PHRED_MAX_STRING_LEN )
    {
      pstrcat( chromatData->comment, "\n", PHRED_MAX_STRING_LEN );
    }
  }

  for( i = 0; i < numEntry; ++i )
  {
    if( ce[i].used == 0 &&
        ( strlen( chromatData->comment ) +
          strlen( ce[i].label ) +
          strlen( ce[i].value ) +
          2 ) < PHRED_MAX_STRING_LEN )
    {
      pstrcat( chromatData->comment, ce[i].label, PHRED_MAX_STRING_LEN );
      pstrcat( chromatData->comment, "=", PHRED_MAX_STRING_LEN );
      pstrcat( chromatData->comment, ce[i].value, PHRED_MAX_STRING_LEN );
      pstrcat( chromatData->comment, "\n", PHRED_MAX_STRING_LEN );
    }
  }

  return( OK );
}

/*
** status:
**
**   0 = OK
**   1 = file reading error
**   2 = no trace (and no bases assumed)
**   3 = no bases (but there is trace)
**  -1 = fatal error
*/

/*
** Read SCF file format >= 3.0.
*/
#ifdef ANSI_C
int readSCF3( char *fn, FILE *fp, SCFHeader *header, ChromatData *chromatData, int numPoint, int numBase, int *status )
#else
int readSCF3( fn, fp, header, chromatData, numPoint, numBase, status )
char *fn;
FILE *fp;
SCFHeader *header;
ChromatData *chromatData;
int numPoint;
int numBase;
int *status;
#endif
{
  int i, j;
  int numByte;
  int numCommEntry;
  uint4 *index;
  char *comments;
  CommentEntry *ce;
  Option *option;

  option = getOption();

  *status = 0;

  /*
  ** Read comments.
  */
  comments = (char *)malloc( header->comments_size + 1 );
  if( comments == NULL )
  {
    fprintf( stderr, "readSCF3: error: unable to allocate memory\n" );
    *status = -1;
    return( ERROR );
  }

  if( header->comments_size > 0 )
  {
    if( fseek( fp, header->comments_offset, 0 ) != 0 )
    {
      fprintf( stderr, "readSCF3: bad status: fseek\n" );
      free( comments );
      *status = 1;
      return( ERROR );
    }
    if( fread( comments, header->comments_size, 1, fp ) == 0 )
    {
      fprintf( stderr, "readSCF3: bad file read\n" );
      free( comments );
      *status = 1;
      return( ERROR );
    }
  }

  /*
  ** Extract comments.
  */
  ce = parseSCFComment( header->comments_size, comments, &numCommEntry );
  if( ce == NULL )
  {
    fprintf( stderr, "readSCF3: error: bad status: parseSCFComment\n" );
    free( comments );
    *status = -1;
    return( ERROR );
  }

  readComment3( numCommEntry, ce, chromatData );

  /*
  ** Free comment memory.
  */
  free( ce );
  free( comments );

  /*
  ** Initialize thumbprint.
  */
  memset( chromatData->thumbPrint, 0, 10 );

  if( numPoint == 0 )
  {
    *status = 2;
    return( ERROR );
  }

  /*
  ** Read sample points.
  */
  if( fseek( fp, header->samples_offset, 0 ) != 0 )
  {
    fprintf( stderr, "readSCF3: bad status: fseek\n" );
    *status = 2;
    return( ERROR );
  }

  if( header->sample_size == 2 )
  {
    uint2 *data;
    uint2 prev;

    /*
    ** Allocate memory.
    */
    numByte = numPoint * sizeof( uint2 );
    data = (uint2 *)malloc( numByte );
    if( data == NULL )
    {
      fprintf( stderr, "readSCF3: error: unable to allocate memory\n" );
      *status = -1;
      return( ERROR );
    }

    /*
    ** Read traces.
    */
    for( j = 0; j < 4; ++j )
    {
      if( fread( (char *)data, 2, numPoint, fp ) != numPoint )
      {
        fprintf( stderr, "readSCF3: error: unable to read trace\n" );
        free( (char *)data );
        *status = 2;
        return( ERROR );
      }
      prev = 0;
      for( i = 0; i < numPoint; ++i )
      {
        data[i] = inSwpUint2( (char *) &(data[i]) ) + prev;
        prev = data[i];
      }
      prev = 0;
      for( i = 0; i < numPoint; ++i )
      {
        data[i] = data[i] + prev;
        chromatData->trace[j][i] = data[i];
        prev = data[i];
      }
    }

    /*
    ** Find minimum and maximum trace values.
    */
    findTraceExtrema( chromatData );

    free( (char *)data );
  }
  else
  {
    uint1 *data;
    uint1 prev;

    /*
    ** Allocate memory.
    */
    numByte = numPoint * sizeof( uint1 );
    data = (uint1 *)malloc( numByte );
    if( data == NULL )
    {
      fprintf( stderr, "readSCF3: error: unable to allocate memory\n" );
      *status = -1;
      return( ERROR );
    }

    /*
    ** Read traces.
    */
    for( j = 0; j < 4; ++j )
    {
      if( fread( (char *)data, 1, numPoint, fp ) != numPoint )
      {
        fprintf( stderr, "readSCF3: error: unable to read trace\n" );
        free( (char *)data );
        *status = 2;
        return( ERROR );
      }
      prev = 0;
      for( i = 0; i < numPoint; ++i )
      {
        data[i] = data[i] + prev;
        prev = data[i];
      }
      prev = 0;
      for( i = 0; i < numPoint; ++i )
      {
        data[i] = data[i] + prev;
        chromatData->trace[j][i] = data[i];
        prev = data[i];
      }
    }

    /*
    ** Find minimum and maximum trace values.
    */
    findTraceExtrema( chromatData );

    free( (char *)data );
  }

  /*
  ** Read bases.
  */
  if( fseek( fp, header->bases_offset + numBase * 8, 0 ) != 0 )
  {
    fprintf( stderr, "readSCF3: bad status: fseek\n" );
    *status = 3;
    return( ERROR );
  }
  if( fread( chromatData->base, 1, numBase, fp ) != numBase )
  {
    fprintf( stderr, "readSCF3: error: unable to read bases\n" );
    *status = 3;
    return( ERROR );
  }

  /*
  ** Read base indices.
  */
  numByte = ( numBase + 1 ) * sizeof( uint4 );
  index = (uint4 *)malloc( numByte );
  if( index == NULL )
  {
    fprintf( stderr, "readSCF3: error: unable to allocate memory\n" );
    *status = -1;
    return( ERROR );
  }
  if( fseek( fp, header->bases_offset, 0 ) != 0 )
  {
    fprintf( stderr, "readSCF3: bad status: fseek\n" );
    free( (char *)index );
    *status = 3;
    return( ERROR );
  }
  if( fread( (char *)index, 4, numBase, fp ) != numBase )
  {
    fprintf( stderr, "readSCF3: error: unable to read base indices\n" );
    free( (char *)index );
    *status = 3;
    return( ERROR );
  }
  for( i = 0; i < numBase; ++i )
  {
    chromatData->baseLoc[i]  = inSwpUint4( (char *) &(index[i]) );
    chromatData->baseQual[i] = 0;
  }
  free( (char *)index );

  if( option->verboseOption && option->verboseLevel >= 16 )
  {
    fprintf( stderr,
             "readSCF3: machine name:     %s\n",
             chromatData->machineName );
    fprintf( stderr,
             "readSCF3: gel name:         %s\n",
             chromatData->gelName );
    fprintf( stderr,
             "readSCF3: sample name:      %s\n",
             chromatData->sampleName );
    fprintf( stderr,
             "readSCF3: primer ID:        %s\n",
             chromatData->primerID );
    fprintf( stderr,
             "readSCF3: lane number:      %d\n",
             chromatData->laneNumber );
    fprintf( stderr,
             "readSCF3: comment:          %s\n",
             chromatData->comment );
    fprintf( stderr,
             "readSCF3: number of scans:  %d\n",
             chromatData->numPoint );
    fprintf( stderr,
             "readSCF3: sample size:      %d bytes\n",
             header->sample_size );
  }

  *status = 0;

  return( OK );
}

