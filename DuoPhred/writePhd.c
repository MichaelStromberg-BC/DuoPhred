/** writePhd.c **/

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
int writePhd( char *filename, char *seqName, Option *option, PhredData *phredData )
#else
int writePhd( filename, seqName, option, phredData )
char *filename;
char *seqName;
Option *option;
PhredData *phredData;
#endif
{
  int i;
  int set;
  int numBase;
  int pos;
  int qual;
  char base;
  char base2;
  int qual2;
  int pos2;
  FILE *fp;

  set = phredData->dataSet[0];
  numBase = phredData->numBase[set];

  /*
  ** Open output file.
  */
  if( ( fp = fopen( filename, "w+" ) ) == NULL )
  {
    fprintf( stderr, "unable to open file %s\n", filename );
    return( ERROR );
  }

  /*
  ** Write sequence delimiter - first line.
  */
  fprintf( fp, "BEGIN_SEQUENCE %s\n", seqName );
  fprintf( fp, "\n" );

  /*
  ** Write comment section delimiters.
  */
  fprintf( fp, "BEGIN_COMMENT\n" );

  fprintf( fp, "\n" );

  fprintf( fp, "CHROMAT_FILE: %s\n", getFileName( option->readFileName ) );

  if( ( ( phredData->thumbPrint[0] & 0xff ) != 0 ) ||
      ( ( phredData->thumbPrint[1] & 0xff ) != 0 ) ||
      ( ( phredData->thumbPrint[2] & 0xff ) != 0 ) ||
      ( ( phredData->thumbPrint[3] & 0xff ) != 0 ) ||
      ( ( phredData->thumbPrint[4] & 0xff ) != 0 ) ||
      ( ( phredData->thumbPrint[5] & 0xff ) != 0 ) ||
      ( ( phredData->thumbPrint[6] & 0xff ) != 0 ) ||
      ( ( phredData->thumbPrint[7] & 0xff ) != 0 ) ||
      ( ( phredData->thumbPrint[8] & 0xff ) != 0 ) ||
      ( ( phredData->thumbPrint[9] & 0xff ) != 0 ) )
  {
    fprintf( fp, "ABI_THUMBPRINT: %3.3o%3.3o%3.3o%3.3o%3.3o%3.3o%3.3o%3.3o%3.3o%3.3o\n",
                 phredData->thumbPrint[0] & 0xff,
                 phredData->thumbPrint[1] & 0xff,
                 phredData->thumbPrint[2] & 0xff,
                 phredData->thumbPrint[3] & 0xff,
                 phredData->thumbPrint[4] & 0xff,
                 phredData->thumbPrint[5] & 0xff,
                 phredData->thumbPrint[6] & 0xff,
                 phredData->thumbPrint[7] & 0xff,
                 phredData->thumbPrint[8] & 0xff,
                 phredData->thumbPrint[9] & 0xff );
  }
  else
  {
    fprintf( fp, "ABI_THUMBPRINT: 0\n" );
  }

  fprintf( fp, "PHRED_VERSION: %s\n", getVersion() );
  fprintf( fp, "CALL_METHOD: %s\n", option->call ? "phred" : "abi" );
  fprintf( fp, "QUALITY_LEVELS: 99\n" );
  fprintf( fp, "TIME: %s\n", getTime() );
  fprintf( fp, "TRACE_ARRAY_MIN_INDEX: 0\n" );
  fprintf( fp, "TRACE_ARRAY_MAX_INDEX: %d\n", phredData->numPoint - 1 );

  if( option->trimPHDData == 0 )
  {
    fprintf( fp,
             "TRIM: %d %d %5.4f\n",
             phredData->lftPhredTrim,
             phredData->rhtPhredTrim,
             phredData->trimSetValue );
  }
  else
  {
    fprintf( fp,
             "TRIM: %d %d -1.00\n",
             0,
             phredData->numBase[set] -
             phredData->rghtTrimPoint -
             phredData->leftTrimPoint - 1 );
  }

  /*
  ** Chemistry type.
  */
  if( phredData->chemType == PRIMER_CHEM )
  {
    fprintf( fp, "CHEM: prim\n" );
  }
  else
  if( phredData->chemType == TERMINATOR_CHEM )
  {
    fprintf( fp, "CHEM: term\n" );
  }
  else
  if( phredData->chemType == UNKNOWN_CHEM )
  {
    fprintf( fp, "CHEM: unknown\n" );
  }
  else
  {
    fprintf( stderr, "writePhd: error: unknow chemistry code: %d\n", phredData->chemType );
    fprintf( fp, "CHEM: unknown\n" );
  }

  /*
  ** Dye type.
  */
  if( phredData->dyeType == RHODAMINE_DYE )
  {
    fprintf( fp, "DYE: rhod\n" );
  }
  else
  if( phredData->dyeType == BIG_DYE_DYE )
  {
    fprintf( fp, "DYE: big\n" );
  }
  else
  if( phredData->dyeType == ENERGY_TRANSFER_DYE )
  {
    fprintf( fp, "DYE: ET\n" );
  }
  else
  if( phredData->dyeType == D_RHODAMINE_DYE )
  {
    fprintf( fp, "DYE: d-rhod\n" );
  }
  else
  if( phredData->dyeType == BODIPY_DYE )
  {
    fprintf( fp, "DYE: bodipy\n" );
  }
  else
  if( phredData->dyeType == UNKNOWN_DYE )
  {
    fprintf( fp, "DYE: unknown\n" );
  }
  else
  {
    fprintf( stderr, "writePhd: error: unknown dye code: %d\n", phredData->dyeType );
    fprintf( fp, "DYE: unknown\n" );
  }

  fprintf( fp, "\n" );

  fprintf( fp, "END_COMMENT\n" );
  fprintf( fp, "\n" );

  /*
  ** Write "DNA" opening section delimiter.
  */
  fprintf( fp, "BEGIN_DNA\n" );

  /*
  ** Write left-cutoff bases
  ** (This will need modification if
  **  editing is re-enabled.)
  */
  if( option->trimPHDData == 0 )
  {
    for( i = 0; i < phredData->leftTrimPoint; ++i )
    {
      base = tolower( phredData->base[set][i] );
      qual = phredData->baseQual[set][i];
      pos  = phredData->baseLoc[set][i];
	  fprintf( fp, "%c %1d %-d\n", base, qual, pos);
      //fprintf( fp, "%c %1d %-d - %c %1d %-d\n", base, qual, pos, base2, qual2, pos2 );
    }
  }

  /*
  ** Write non-cutoff bases.
  ** (This will need modification if
  **  editing is re-enabled.)
  */
  for( i = phredData->leftTrimPoint;
       i < numBase - phredData->rghtTrimPoint;
       ++i )
  {
    base = tolower( phredData->base[set][i] );
    qual = phredData->baseQual[set][i];
    pos  = phredData->baseLoc[set][i];
    fprintf( fp, "%c %1d %-d\n", base, qual, pos );
  }

  /*
  ** Write right-cutoff bases.
  ** (This will need modification if
  **  editing is re-enabled.)
  */
  if( option->trimPHDData == 0 )
  {
    for( i = numBase - phredData->rghtTrimPoint;
         i < numBase;
         ++i )
    {
      base = tolower( phredData->base[set][i] );
      qual = phredData->baseQual[set][i];
      pos  = phredData->baseLoc[set][i];
      fprintf( fp, "%c %1d %-d\n", base, qual, pos );
    }
  }

  /*
  ** Write "DNA" closing section delimiter.
  */
  fprintf( fp, "END_DNA\n" );
  fprintf( fp, "\n" );

  /*
  ** Write file closing delimiter.
  */
  fprintf( fp, "END_SEQUENCE" );

  /*
  ** Write a carriage return.  And test for error.
  ** (Remove phd file on error.)
  */
  if( fwrite( "\n", 1, 1, fp ) != 1 )
  {
    fprintf( stderr,
             "writePhd: error: unable to write to %s\n",
             filename );
    fclose( fp );
    delFile( filename );
    return( ERROR );
  }

  /*
  ** Close PHD file.
  */
  fclose( fp );

  return( OK );
}
