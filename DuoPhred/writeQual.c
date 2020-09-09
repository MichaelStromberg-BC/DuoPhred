/** writeQual.c **/

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
**    * writeQual.c                                                    *     **
**    * benefits from ideas in code written by LaDeana Hillier and     *     **
**    * Tim Gleeson.                                                   *     **
**                                                                           **
*******************************************************************************
*/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "phred.h"

#define CharPerLine 50 /* For output formatting */

#ifdef ANSI_C
static char *tail( char *pathname )
#else
static char *tail( pathname )
char *pathname;
#endif
{
  char *a;
  if ((a = (char *) strrchr(pathname,PATHSEP))==NULL)
  {
    a = pathname;
  }
  else
  {
    a++;
  }
  return( a );
}

/*
** writeScore
**
**
** Write quality values in either FASTA (with cut-off base values) or
** XBAP (without cut-off base values) form.
**
*/
#ifdef ANSI_C
static int writeScore( PhredData *phredData, Option *option, FILE *fp )
#else
static int writeScore( phredData, option, fp )
PhredData *phredData;
Option *option;
FILE *fp;
#endif
{
  int i;
  int set;
  int qual;
  int numBase;
  int lineLen;
  int zero;
  char string[64];

  zero = 0;

  lineLen = 0;

  set = phredData->dataSet[0];
  numBase = phredData->numBase[set];

  if( option->qualType == 0 ||
      option->qualType == 1 )
  {
    if( option->qualType == 0 && option->trimFastaData == 0 )
    {
      /*
      ** write left-cutoff as zero quality value assignments
      */
      for( i = 0;
           i < phredData->leftTrimPoint;
           ++i )
      {
        sprintf( string, "%1d ", zero );
        fprintf( fp, "%s", string );
        lineLen += strlen( string );
        if( lineLen >= CharPerLine )
        {
          if( fputc( '\n', fp ) == EOF && ferror( fp ) != 0 )
          {
            fprintf( stderr,
                     "writeScore: unable to write to qual file\n" );
            return( ERROR );
          }
          lineLen = 0;
        }
      }
    }

    /*
    ** write non-cutoff calls
    */
    for( i = phredData->leftTrimPoint;
         i < numBase - phredData->rghtTrimPoint;
         ++i )
    {
      qual = phredData->baseQual[set][i];
      sprintf( string, "%1d ", qual );
      fprintf( fp, "%s", string );
      lineLen += strlen( string );
      if( lineLen >= CharPerLine )
      {
        if( fputc( '\n', fp ) == EOF && ferror( fp ) != 0 )
        {
          fprintf( stderr,
                   "writeScore: unable to write to qual file\n" );
          return( ERROR );
        }
        lineLen = 0;
      }
    }

    if( option->qualType == 0 && option->trimFastaData == 0 )
    {
      /*
      ** write right-cutoff as zero quality value assignments
      */
      for( i = numBase - phredData->rghtTrimPoint;
           i < numBase;
           ++i )
      {
        sprintf( string, "%1d ", zero );
        fprintf( fp, "%s", string );
        lineLen += strlen( string );
        if( lineLen >= CharPerLine )
        {
          if( fputc( '\n', fp ) == EOF && ferror( fp ) != 0 )
          {
            fprintf( stderr,
                     "writeScore: unable to write to qual file\n" );
            return( ERROR );
          }
          lineLen = 0;
        }
      }
    }

    if( lineLen != 0 )
    {
      if( fputc( '\n', fp ) == EOF && ferror( fp ) != 0 )
      {
        fprintf( stderr,
                 "writeScore: unable to write to qual file\n" );
        return( ERROR );
      }
    }

  }
  else
  if( option->qualType == 2 )
  {
    /*
    ** Write all calls.
    */
    for( i = 0; i < numBase; ++i )
    {
      qual = phredData->baseQual[set][i];
      sprintf( string, "%1d ", qual );
      fprintf( fp, "%s", string );
      lineLen += strlen( string );
      if( lineLen >= CharPerLine )
      {
        if( fputc( '\n', fp ) == EOF && ferror( fp ) != 0 )
        {
          fprintf( stderr,
                   "writeScore: unable to write to qual file\n" );
          return( ERROR );
        }
        lineLen = 0;
      }
    }
    if( lineLen != 0 )
    {
      if( fputc( '\n', fp ) == EOF && ferror( fp ) != 0 )
      {
        fprintf( stderr,
                 "writeScore: unable to write to qual file\n" );
        return( ERROR );
      }
    }
  }

  return( OK );

} /* writeScore */



/*
** writeQual
**
** Write quality file.
**
*/
#ifdef ANSI_C
int writeQual( char *filename, char *seqName, Option *option, PhredData *phredData )
#else
int writeQual( filename, seqName, option, phredData )
char *filename;
char *seqName;
Option *option;
PhredData *phredData;
#endif
{
  int set;
  int numBase;
  int type;
  FILE *fp;

#ifdef OMIT_EMPTY_ENTRY
  /*
  ** Do not write sequences with no bases.
  */
  set = phredData->dataSet[0];
  numBase = phredData->numBase[set];

  if( numBase == 0 )
  {
    return( OK );
  }
#endif

  if( option->qaOption == 0 )
  {
    /*
    ** Open quality file.
    */
    if( ( fp = fopen( filename, "w+" ) ) == NULL )
    {
      fprintf( stderr, "writeData: unable to open file %s\n", filename );
      return( ERROR );
    }
  }
  else
  {
    fp = option->qualFP;
  }

  /*
  ** Write header.
  */
  if( option->qualType == 1 )
  {
    type = 1;
  }
  else
  {
    type = 0;
  }
  if( writeHeader( phredData, seqName, option, fp, type ) == ERROR )
  {
    if( option->qaOption == 0 )
    {
      /*
      ** close quality file
      */
      fclose( fp );
    }
    return( ERROR );
  }

/*
** set scores to divide calls into ranges for
** evaluating PHRED using cross_match
** setRange( phredData );
*/

  /*
  ** write quality file
  */
  if( writeScore( phredData, option, fp ) == ERROR )
  {
    if( option->qaOption == 0 )
    {
      /*
      ** close quality file
      */
      fclose( fp );
    }
    return( ERROR );
  }

  if( option->qaOption == 0 )
  {
    /*
    ** close quality file
    */
    fclose( fp );
  }

  return( 0 );
}



/*
** Adjust quality scores to identify base
** location from left cutoff.
*/
#ifdef ANSI_C
int setRange( PhredData *phredData )
#else
int setRange( phredData )
PhredData *phredData;
#endif
{
  int i;
  int ic;
  int set;
  int numBase;
  int rngAdd;
  int rngwid;

  set = phredData->dataSet[0];
  numBase = phredData->numBase[set];

  rngwid = 100;

  rngAdd = -10;
  ic = 0;
  for( i = phredData->leftTrimPoint; i < numBase; ++i )
  {
    if( ( ic % rngwid ) == 0 && rngAdd < 120 )
    {
      rngAdd += 10;
    }
    phredData->baseQual[set][i] += rngAdd;
    ++ic;
  }

  return( OK );
}


