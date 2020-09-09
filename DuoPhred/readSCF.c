/** readSCF.c **/

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
**    * readSCF.c                                                      *     **
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
#include "typeDef.h"
#include "rwUtil.h"
#include "chromatData.h"
#include "freeChromatData.h"
#include "readSCF.h"


#ifdef ANSI_C
static int readSCFHeader( FILE *fp, SCFHeader *h )
#else
static int readSCFHeader( fp, h )
FILE *fp;
SCFHeader *h;
#endif
{
  if( readUint4( fp, &h->magic_number )     == ERROR ) return( ERROR );
  if( readUint4( fp, &h->samples )          == ERROR ) return( ERROR );
  if( readUint4( fp, &h->samples_offset )   == ERROR ) return( ERROR );
  if( readUint4( fp, &h->bases )            == ERROR ) return( ERROR );
  if( readUint4( fp, &h->bases_left_clip )  == ERROR ) return( ERROR );
  if( readUint4( fp, &h->bases_right_clip ) == ERROR ) return( ERROR );
  if( readUint4( fp, &h->bases_offset )     == ERROR ) return( ERROR );
  if( readUint4( fp, &h->comments_size )    == ERROR ) return( ERROR );
  if( readUint4( fp, &h->comments_offset )  == ERROR ) return( ERROR );
  if( fread( h->version, sizeof( h->version ), 1, fp) != 1 ) return( ERROR );
  if( readUint4( fp, &h->sample_size )      == ERROR ) return( ERROR );
  if( readUint4( fp, &h->code_set )         == ERROR ) return( ERROR );
  return( OK );
}



/*
** Read the SCF format sequence with name `fn' into `seq'.
*/

/*
** status:
**
**   0 = OK
**   1 = file reading error
**   2 = no trace (and no bases assumed)
**   3 = no bases (but there is trace)
**  -1 = fatal error
*/

#ifdef ANSI_C
ChromatData *readSCF( char *fn, int *status )
#else
ChromatData *readSCF( fn, status )
char *fn;
int *status;
#endif
{
  int numBase;
  int numPoint;
  int versionSwitch;
  int lstat;
  SCFHeader header;
  ChromatData *chromatData;
  FILE *fp;

  /*
  ** Open file for reading.
  */
  fp = fopen( fn, "rb" );
  if( fp == NULL )
  {
    fprintf( stderr, "readSCF: unable to open file %s\n", fn );
    *status = 1;
    return( NULL );
  }

  /*
  ** Read header.
  */
  if( readSCFHeader( fp, &header ) == ERROR )
  {
    fprintf( stderr, "readSCF: unable to read %s header\n", fn );
    fclose( fp );

    numPoint = 0;
    numBase  = 0;

    chromatData = allocChromatData( numPoint, numBase );
    if( chromatData == NULL )
    {
      fprintf( stderr, "readSCF: unable to allocate memory\n" );
      fclose( fp );
      *status = -1;
      return( NULL );
    }

    chromatData->fileType           = SCFFormat;
    chromatData->primerLoc          = 0;
    chromatData->avgSpacing         = 0.0;
    chromatData->machineName[0]     = '\0';
    chromatData->sampleName[0]      = '\0';
    chromatData->primerID[0]        = '\0';
    chromatData->signalStrength[0]  = 0;
    chromatData->signalStrength[1]  = 0;
    chromatData->signalStrength[2]  = 0;
    chromatData->signalStrength[3]  = 0;
    chromatData->gelName[0]         = '\0';
    chromatData->laneNumber         = 0;
    chromatData->processing[0]      = '\0';
    chromatData->reTracker[0]       = '\0';
    chromatData->comment[0]         = '\0';
    chromatData->convProg[0]        = '\0';
    chromatData->source[0]          = '\0';
    pstrcpy( chromatData->fileName, fn, PHRED_PATH_MAX );
    memset( chromatData->thumbPrint, 0, 10 * sizeof( char ) );

    *status = 1;
    return( chromatData );
  }

  numPoint = header.samples;
  numBase  = header.bases;

  /*
  ** Allocate chromatData memory.
  */
  chromatData = allocChromatData( numPoint, numBase );
  if( chromatData == NULL )
  {
    fprintf( stderr, "readSCF: unable to allocate memory\n" );
    fclose( fp );
    *status = -1;
    return( NULL );
  }

  /*
  ** Initialize values.
  */
  chromatData->fileType           = SCFFormat;
  chromatData->primerLoc          = 0;
  chromatData->avgSpacing         = 0.0;
  chromatData->machineName[0]     = '\0';
  chromatData->sampleName[0]      = '\0';
  chromatData->primerID[0]        = '\0';
  chromatData->signalStrength[0]  = 0;
  chromatData->signalStrength[1]  = 0;
  chromatData->signalStrength[2]  = 0;
  chromatData->signalStrength[3]  = 0;
  chromatData->gelName[0]         = '\0';
  chromatData->laneNumber         = 0;
  chromatData->processing[0]      = '\0';
  chromatData->reTracker[0]       = '\0';
  chromatData->comment[0]         = '\0';
  chromatData->convProg[0]        = '\0';
  chromatData->source[0]          = '\0';
  pstrcpy( chromatData->fileName, fn, PHRED_PATH_MAX );
  memset( chromatData->thumbPrint, 0, 10 * sizeof( char ) );

  /*
  ** Check the file version.
  */
  if( strtod( header.version, (char **)NULL ) < 2.0 )
  {
    versionSwitch = 2;
    header.sample_size = 1;
  }
  else
  if( strtod( header.version, (char **)NULL ) < 3.0 )
  {
    versionSwitch = 2;
  }
  else
  {
    versionSwitch = 3;
  }

  if( versionSwitch == 2 )
  {
    if( readSCF2( fn, fp, &header, chromatData, numPoint, numBase, &lstat ) == ERROR )
    {
      if( lstat == 1 || lstat == 2 )
      {
        chromatData->numPoint = 0;
        chromatData->numBase  = 0;
        *status = lstat;
      }
      else
      if( lstat == 3 )
      {
        chromatData->numBase = 0;
        *status = lstat;
      }
      else
      if( lstat == 0 )
      {
        fprintf( stderr,
                 "readSCF: internal inconsistency: readSCF2 returns ERROR with status OK\n" );
        chromatData->numPoint = 0;
        chromatData->numBase  = 0;
        *status = -1;
      }
      else
      if( lstat == -1 )
      {
        chromatData->numPoint = 0;
        chromatData->numBase  = 0;
        *status = lstat;
      }
      else
      {
        fprintf( stderr,
                 "readSCF: unknown status: %d\n", lstat );
        chromatData->numPoint = 0;
        chromatData->numBase  = 0;
        *status = -1;
      }

      fclose( fp );
      return( chromatData );
    }
  }
  else
  if( versionSwitch == 3 )
  {
    if( readSCF3( fn, fp, &header, chromatData, numPoint, numBase, &lstat ) == ERROR )
    {
      if( lstat == 1 || lstat == 2 )
      {
        chromatData->numPoint = 0;
        chromatData->numBase  = 0;
        *status = lstat;
      }
      else
      if( lstat == 3 )
      {
        chromatData->numBase = 0;
        *status = lstat;
      }
      else
      if( lstat == 0 )
      {
        fprintf( stderr,
                 "readSCF: internal inconsistency: readSCF2 returns ERROR with status OK\n" );
        chromatData->numPoint = 0;
        chromatData->numBase  = 0;
        *status = -1;
      }
      else
      if( lstat == -1 )
      {
        chromatData->numPoint = 0;
        chromatData->numBase  = 0;
        *status = lstat;
      }
      else
      {
        fprintf( stderr,
                 "readSCF: unknown status: %d\n", lstat );
        chromatData->numPoint = 0;
        chromatData->numBase  = 0;
        *status = -1;
      }

      fclose( fp );
      return( chromatData );
    }
  }

  *status = lstat;

  /*
  ** Finished with the file.
  */
  fclose( fp );

  return( chromatData );
}

