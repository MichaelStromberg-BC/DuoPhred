/** writeSeq.c **/

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
**    * wrietSeq.c                                                     *     **
**    * benefits from ideas in code written by LaDeana Hillier and     *     **
**    * Tim Gleeson.                                                   *     **
**                                                                           **
*******************************************************************************
*/

/*
** Write sequence file.
*/

#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "phred.h"

#define BasesPerLine 50

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
** writeHeader
**
** Write a header in either FASTA or XBAP form.
**
*/
#ifdef ANSI_C
int writeHeader( PhredData *phredData, char *seqName, Option *option, FILE *fp, int type )
#else
int writeHeader( phredData, seqName, option, fp, type )
PhredData *phredData;
char *seqName;
Option *option;
FILE *fp;
int type;
#endif
{
  int set;
  int left;
  int numBase;
  int numTrim;

  set = phredData->dataSet[0];
  numBase = phredData->numBase[set];

  if( type == 0 )
  {
    if( option->trimFastaData == 1 )
    {
      numBase = phredData->numBase[set] - phredData->rghtTrimPoint - phredData->leftTrimPoint;
      numTrim = numBase;
      left = 0;
    }
    else
    {
      numBase = phredData->numBase[set];
      numTrim = phredData->numBase[set] - phredData->rghtTrimPoint - phredData->leftTrimPoint;
      left    = phredData->leftTrimPoint;;
    }

    fprintf( fp, ">%s %6d %6d %6d %4s",
                 tail( seqName ),
                 numBase,
                 left,
                 numTrim,
                 ( phredData->fileType == ABIFormat ) ? "ABI":
                 ( phredData->fileType == SCFFormat ) ? "SCF":
                 ( phredData->fileType == MD1Format ) ? "ESD":
                 ( phredData->fileType == MD2Format ) ? "ESD":
                 "   "
           );

    if( option->trimFastaData == 1 )
    {
      fprintf( fp, " trimmed\n" );
    }
    else
    {
      fprintf( fp, "\n" );
    }
  }
  else
  {
    fprintf( fp, ";%6d%6d%6d%-4s%-12s\n",
                 numBase,
                 phredData->leftTrimPoint,
                 numBase - phredData->rghtTrimPoint - phredData->leftTrimPoint,
                 ( phredData->fileType == ABIFormat ) ? "ABI":
                 ( phredData->fileType == SCFFormat ) ? "SCF":
                 ( phredData->fileType == MD1Format ) ? "ESD":
                 ( phredData->fileType == MD2Format ) ? "ESD":
                 "   ",
                 tail( seqName ) );
  }

  return( OK );

} /* write_header */


/*
** writeCutoff
**
** Write cut-off bases as XBAP comments.
**
*/
#ifdef ANSI_C
static int writeCutoff( PhredData *phredData, FILE *fp )
#else
static int writeCutoff( phredData, fp )
PhredData *phredData;
FILE *fp;
#endif
{
  int i;
  int lineLen;
  int numBase;
  int set;
  char base;

  set = phredData->dataSet[0];
  numBase = phredData->numBase[set];

  /*
  ** output left cut off
  */
  lineLen = 0;
  for( i = 0; i < phredData->leftTrimPoint; i++ )
  {
    if( !lineLen )
    {
       fprintf( fp, ";<" );
    }
    base = phredData->base[set][i];
    if( fputc( base, fp ) == EOF && ferror( fp ) != 0 )
    {
      fprintf( stderr,
               "writeCutoff: unable to write to seq file\n" );
      return( ERROR );
    }
    if( ++lineLen == BasesPerLine )
    {
      if( fputc( '\n', fp ) == EOF && ferror( fp ) != 0 )
      {
        fprintf( stderr,
                 "writeCutoff: unable to write to seq file\n" );
        return( ERROR );
      }
      lineLen = 0;
    }
  }
  if( lineLen )
  {
    if( fputc( '\n', fp ) == EOF && ferror( fp ) != 0 )
    {
      fprintf( stderr,
               "writeCutoff: unable to write to seq file\n" );
      return( ERROR );
    }
  }

  /*
  ** output right cut off
  */
  lineLen = 0;
  for( i = numBase - phredData->rghtTrimPoint; i < numBase; i++ )
  {
    if( !lineLen )
    {
      fprintf( fp, ";>" );
    }
    base = phredData->base[set][i];
    if( fputc( base, fp ) == EOF && ferror( fp ) != 0 )
    {
      fprintf( stderr,
               "writeCutoff: unable to write to seq file\n" );
      return( ERROR );
    }
    if( ++lineLen == BasesPerLine )
    {
      if( fputc( '\n', fp ) == EOF && ferror( fp ) != 0 )
      {
        fprintf( stderr,
                 "writeCutoff: unable to write to seq file\n" );
        return( ERROR );
      }
      lineLen = 0;
    }
  }
  if( lineLen )
  {
    if( fputc( '\n', fp ) == EOF && ferror( fp ) != 0 )
    {
      fprintf( stderr,
               "writeCutoff: unable to write to seq file\n" );
      return( ERROR );
    }
  }

  return( OK );
}

/*
** writeBase
**
** Write bases in either FASTA (with cut-off bases) or
** XBAP (without cut-off bases) form.
*/
#ifdef ANSI_C
static int writeBase( PhredData *phredData, Option *option, FILE *fp )
#else
static int writeBase( phredData, option, fp )
PhredData *phredData;
Option *option;
FILE *fp;
#endif
{
  int i;
  int set;
  int numBase;
  int lineLen;
  char base;
  lineLen = 0;

  set = phredData->dataSet[0];
  numBase = phredData->numBase[set];

  if( option->seqType == 0 && option->trimFastaData == 0 )
  {
    /*
    ** Write left-cutoff as low quality calls.
    */
    for( i = 0; i < phredData->leftTrimPoint; ++i )
    {
      base = phredData->base[set][i];
      if( fputc( base, fp ) == EOF && ferror( fp ) != 0 )
      {
        fprintf( stderr,
                 "writeCutoff: unable to write to seq file\n" );
        return( ERROR );
      }
      ++lineLen;
      if( lineLen == BasesPerLine )
      {
        if( fputc( '\n', fp ) == EOF && ferror( fp ) != 0 )
        {
          fprintf( stderr,
                   "writeCutoff: unable to write to seq file\n" );
          return( ERROR );
        }
        lineLen = 0;
      }
    }
  }

  /*
  ** Write non-cutoff calls.
  */
  for( i = phredData->leftTrimPoint; i < numBase - phredData->rghtTrimPoint; ++i )
  {
    base = phredData->base[set][i];
    if( fputc( base, fp ) == EOF && ferror( fp ) != 0 )
    {
      fprintf( stderr,
               "writeCutoff: unable to write to seq file\n" );
      return( ERROR );
    }
    if( ++lineLen == BasesPerLine )
    {
      if( fputc( '\n', fp ) == EOF && ferror( fp ) != 0 )
      {
        fprintf( stderr,
                 "writeCutoff: unable to write to seq file\n" );
        return( ERROR );
      }
      lineLen = 0;
    }
  }

  if( option->seqType == 0 && option->trimFastaData == 0 )
  {
    /*
    ** Write right-cutoff low quality bases.
    */
    for( i = numBase - phredData->rghtTrimPoint; i < numBase; ++i )
    {
      base = phredData->base[set][i];
      if( fputc( base, fp ) == EOF && ferror( fp ) != 0 )
      {
        fprintf( stderr,
                 "writeCutoff: unable to write to seq file\n" );
        return( ERROR );
      }
      ++lineLen;
      if( lineLen == BasesPerLine )
      {
        if( fputc( '\n', fp ) == EOF && ferror( fp ) != 0 )
        {
          fprintf( stderr,
                   "writeCutoff: unable to write to seq file\n" );
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
               "writeCutoff: unable to write to seq file\n" );
      return( ERROR );
    }
  }

  return( OK );

} /* writeBase */



/*
** writeSeq
**
** Write sequence file.
**
*/
#ifdef ANSI_C
int writeSeq( char *filename, char *seqName, Option *option, PhredData *phredData )
#else
int writeSeq( filename, seqName, option, phredData )
char *filename;
char *seqName;
Option *option;
PhredData *phredData;
#endif
{
  int set;
  int numBase;
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

  if( option->saOption == 0 )
  {
    /*
    ** Open sequence file.
    */
    if( ( fp = fopen( filename, "w+" ) ) == NULL )
    {
      fprintf( stderr, "writeData: unable to open file %s\n", filename );
      return( ERROR );
    }
  }
  else
  {
    fp = option->seqFP;
  }

  /*
  ** Write header.
  */
  if( writeHeader( phredData, seqName, option, fp, option->seqType ) == ERROR )
  {
    if( option->saOption == 0 )
    {
      fclose( fp );
    }
    return( ERROR );
  }

  /*
  ** Write trimmed off bases for xbap.
  */
  if( option->seqType == 1 )
  {
    if( writeCutoff( phredData, fp ) == ERROR )
    {
      if( option->saOption == 0 )
      {
        fclose( fp );
      }
      return( ERROR );
    }
  }

  /*
  ** Write sequence.
  */
  if( writeBase( phredData, option, fp ) == ERROR )
  {
    if( option->saOption == 0 )
    {
      fclose( fp );
    }
    return( ERROR );
  }

  if( option->saOption == 0 )
  {
    /*
    ** Close sequence file.
    */
    fclose( fp );
  }

  return( 0 );
}

