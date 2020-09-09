/** readSCF2.c **/

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
**    * readSCF2.c                                                      *     **
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
static int readSCFSample1( FILE *fp, Samples1 *s )
#else
static int readSCFSample1( fp, s )
FILE *fp;
Samples1 *s;
#endif
{
  if( readUint1( fp, &s->sample_A ) == ERROR ) return( ERROR );
  if( readUint1( fp, &s->sample_C ) == ERROR ) return( ERROR );
  if( readUint1( fp, &s->sample_G ) == ERROR ) return( ERROR );
  if( readUint1( fp, &s->sample_T ) == ERROR ) return( ERROR );
  return( OK );
}


#ifdef ANSI_C
static int readSCFSample2( FILE *fp, Samples2 *s )
#else
static int readSCFSample2( fp, s )
FILE *fp;
Samples2 *s;
#endif
{
  if( readUint2( fp, &s->sample_A ) == ERROR ) return( ERROR );
  if( readUint2( fp, &s->sample_C ) == ERROR ) return( ERROR );
  if( readUint2( fp, &s->sample_G ) == ERROR ) return( ERROR );
  if( readUint2( fp, &s->sample_T ) == ERROR ) return( ERROR );
  return( OK );
}


#ifdef ANSI_C
static int readSCFBase( FILE *fp, Bases *b )
#else
static int readSCFBase( fp, b )
FILE *fp;
Bases *b;
#endif
{
  if( readUint4( fp, &b->peak_index )    == ERROR ) return( ERROR );
  if( readUint1( fp, &b->prob_A )        == ERROR ) return( ERROR );
  if( readUint1( fp, &b->prob_C )        == ERROR ) return( ERROR );
  if( readUint1( fp, &b->prob_G )        == ERROR ) return( ERROR );
  if( readUint1( fp, &b->prob_T )        == ERROR ) return( ERROR );
  if( readUint1( fp, (uint1 *)&b->base ) == ERROR ) return( ERROR );
  if( readUint1( fp, &b->spare[0] )      == ERROR ) return( ERROR );
  if( readUint1( fp, &b->spare[1] )      == ERROR ) return( ERROR );
  if( readUint1( fp, &b->spare[2] )      == ERROR ) return( ERROR );
  return( OK );
}


#ifdef ANSI_C
static int readComment2( int numEntry, CommentEntry *ce, ChromatData *chromatData )
#else
static int readComment2( numEntry, ce, chromatData )
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
** Read SCF file format < 3.0.
*/
#ifdef ANSI_C
int readSCF2( char *fn, FILE *fp, SCFHeader *header, ChromatData *chromatData, int numPoint, int numBase, int *status )
#else
int readSCF2( fn, fp, header, chromatData, numPoint, numBase, status )
char *fn;
FILE *fp;
SCFHeader *header;
ChromatData *chromatData;
int numPoint;
int numBase;
int *status;
#endif
{
  int i;
  int numCommEntry;
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
    fprintf( stderr, "readSCF2: error: unable to allocate memory\n" );
    *status = -1;
    return( ERROR );
  }

  if( header->comments_size > 0 )
  {
    if( fseek( fp, header->comments_offset, 0 ) != 0 )
    {
      fprintf( stderr, "readSCF2: bad status: fseek\n" );
      free( comments );
      *status = 1;
      return( ERROR );
    }
    if( fread( comments, header->comments_size, 1, fp ) == 0 )
    {
      fprintf( stderr, "readSCF2: bad file read\n" );
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
    fprintf( stderr, "readSCF2: error: bad status: parseSCFComment\n" );
    free( comments );
    *status = -1;
    return( ERROR );
  }

  readComment2( numCommEntry, ce, chromatData );

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
    fprintf( stderr, "readSCF2: bad status: fseek\n" );
    chromatData->numPoint = 0;
    chromatData->numBase  = 0;
    *status = 2;
    return( ERROR );
  }
  if( header->sample_size == 2 )
  {
    Samples2 sample;

    for( i = 0; i < numPoint; i++ )
    {
      if( readSCFSample2( fp, &sample ) == ERROR )
      {
        fprintf( stderr, "readSCF2: unable to read sample2 size\n" );
        *status = 2;
        return( ERROR );
      }

      chromatData->trace[0][i] = (FLOAT)sample.sample_A;
      chromatData->trace[1][i] = (FLOAT)sample.sample_C;
      chromatData->trace[2][i] = (FLOAT)sample.sample_G;
      chromatData->trace[3][i] = (FLOAT)sample.sample_T;

    } /* for i */

    /*
    ** Find minimum and maximum trace values.
    */
    findTraceExtrema( chromatData );
  }
  else
  {
    Samples1 sample;

    for( i = 0; i < numPoint; i++ )
    {
      if( readSCFSample1( fp, &sample ) == ERROR )
      {
        fprintf( stderr, "readSCF2: unable to read sample1 size\n" );
        chromatData->numPoint = 0;
        chromatData->numBase  = 0;
        *status = 2;
        return( ERROR );
      }

      chromatData->trace[0][i] = (FLOAT)sample.sample_A;
      chromatData->trace[1][i] = (FLOAT)sample.sample_C;
      chromatData->trace[2][i] = (FLOAT)sample.sample_G;
      chromatData->trace[3][i] = (FLOAT)sample.sample_T;

    } /* for i */

    /*
    ** Find minimum and maximum trace values.
    */
    findTraceExtrema( chromatData );
  }

  /*
  ** Read bases
  */
  if( fseek( fp, header->bases_offset, 0 ) != 0 )
  {
    fprintf( stderr, "readSCF2: bad status: fseek\n" );
    chromatData->numBase = 0;
    *status = 3;
    return( ERROR );
  }
  for( i = 0; i < numBase; i++ )
  {
    Bases base;
    if( readSCFBase( fp, &base ) == ERROR )
    {
      fprintf( stderr, "readSCF2: unable to read base\n" );
      chromatData->numBase = 0;
      *status = 3;
      return( ERROR );
    }
    chromatData->base[i]     = (int)base.base;
    chromatData->baseLoc[i]  = (int)base.peak_index;
    chromatData->baseQual[i] = 0;
  }

  if( option->verboseOption && option->verboseLevel >= 16 )
  {
    fprintf( stderr,
             "readSCF2: machine name:     %s\n",
             chromatData->machineName );
    fprintf( stderr,
             "readSCF2: gel name:         %s\n",
             chromatData->gelName );
    fprintf( stderr,
             "readSCF2: sample name:      %s\n",
             chromatData->sampleName );
    fprintf( stderr,
             "readSCF2: primer ID:        %s\n",
             chromatData->primerID );
    fprintf( stderr,
             "readSCF2: lane number:      %d\n",
             chromatData->laneNumber );
    fprintf( stderr,
             "readSCF2: comment:          %s\n",
             chromatData->comment );
    fprintf( stderr,
             "readSCF2: number of scans:  %d\n",
             chromatData->numPoint );
    fprintf( stderr,
             "readSCF2: sample size:      %d bytes\n",
             header->sample_size );
  }


  *status = 0;

  return( OK );
}

