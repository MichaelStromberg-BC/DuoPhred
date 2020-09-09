/** writeData.c **/

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
** Note: try to have phred write a 'empty' output files
**       on error.
*/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include "phred.h"

extern int errno;

/*
** writeData
**
** Write sequence, quality, and SCF files.
**
*/
#ifdef ANSI_C
int writeData( PhredData *phredData, char *outFileName, Option *option )
#else
int writeData( phredData, outFileName, option )
PhredData *phredData;
char *outFileName;
Option *option;
#endif
{
  int prec;
  char *fileRootName;
  char string[PHRED_PATH_MAX];
  char filename[3*PHRED_PATH_MAX];
  char *seqName;

  /*
  ** Set sequence cutoffs to zero if no trimming.
  */
  if( option->trim == 0 )
  {
    phredData->leftTrimPoint = 0;
    phredData->rghtTrimPoint = 0;
  }
  
  /*
  ** Strip away path specification, if present.
  */
  if( ( fileRootName =
        strrchr( outFileName, PATHSEP ) ) != NULL )
  {
    fileRootName += 1;
  }
  else
  {
    fileRootName = outFileName;
  }

  /*
  ** Sequence name (filename if sequence name is not explicit).
  */
  seqName = option->seqName ? option->seqName : fileRootName;

  /*
  ** Write sequence file, if requested.
  */
  if( option->writeSeq )
  {
    /*
    ** Construct sequence filename.
    */
    if( option->writeSeqName == NULL )
    {
      sprintf( filename,
               "%.*s.seq",
               PHRED_PATH_MAX,
               fileRootName );
    }
    else
    {
      sprintf( filename,
               "%.*s",
               PHRED_PATH_MAX,
               option->writeSeqName );
    }

    /*
    ** Add directory prefix, if requested.
    */
    if( option->sdOption == 1 )
    {
      /*
      ** Linux fails to write filename into filename with sprintf
      */
      pstrcpy( string, filename, PHRED_PATH_MAX );
      sprintf( filename,
               "%.*s%c%.*s",
               PHRED_PATH_MAX,
               option->seqDirName,
               PATHSEP,
               PHRED_PATH_MAX,
               string );
    }

    /*
    ** write sequence file
    */
    if( writeSeq( filename, seqName, option, phredData ) == ERROR )
    {
      return( ERROR );
    }

  } /* if option->writeSeq */

  /*
  ** Write quality file, if requested.
  */
  if( option->writeQual )
  {

    /*
    ** Construct quality filename.
    */
    if( option->writeQualName == NULL )
    {
      sprintf( filename,
               "%.*s.qual",
               PHRED_PATH_MAX,
               fileRootName );
    }
    else
    {
      sprintf( filename,
               "%.*s",
               PHRED_PATH_MAX,
               option->writeQualName );
    }

    /*
    ** Add directory prefix, if requested.
    */
    if( option->qdOption == 1 )
    {
      /*
      ** Linux fails to write filename into filename with sprintf
      */
      pstrcpy( string, filename, PHRED_PATH_MAX );
      sprintf( filename,
               "%.*s%c%.*s",
               PHRED_PATH_MAX,
               option->qualDirName,
               PATHSEP,
               PHRED_PATH_MAX,
               string );
    }

    /*
    ** write quality file
    */
    if( writeQual( filename, seqName, option, phredData ) == ERROR )
    {
      return( ERROR );
    }

  } /* if option->writeQual */

  /*
  ** Write SCF file, if requested.
  */
  if( option->writeScf )
  {
    /*
    ** Construct SCF filename.
    */
    if( option->writeScfName == NULL )
    {
      sprintf( filename,
               "%.*s",
               PHRED_PATH_MAX,
               fileRootName );
    }
    else
    {
      sprintf( filename,
               "%.*s",
               PHRED_PATH_MAX,
               option->writeScfName );
    }
    
    /*
    ** Add directory prefix, if requested.
    */
    if( option->scfDir == 1 )
    {
      /*
      ** Linux fails to write filename into filename with sprintf
      */
      pstrcpy( string, filename, PHRED_PATH_MAX );
      sprintf( filename,
               "%.*s%c%.*s",
               PHRED_PATH_MAX,
               option->scfDirName,
               PATHSEP,
               PHRED_PATH_MAX,
               string );
    }
 
    /*
    ** Check whether SCF file will overwrite input file
    */
    if( compareFiles( filename, option ) == 0 )
    {
      /*
      ** write SCF file
      */
      if( option->scfPrecOption == 1 )
      {
        prec = option->scfPrecision;
      }
      else
      {
        prec = ( (int)phredData->maxTraceValue < 256 ) ? 1 : 2;
      }
      if( option->scfVersion == 3 )
      {
        if( writeSCF3( filename, prec, phredData ) == ERROR )
        {
          return( ERROR );
        }
      }
      else
      {
        if( writeSCF2( filename, prec, phredData ) == ERROR )
        {
          return( ERROR );
        }
      }
    }
    else
    {
      fprintf( stderr, "SCF file not written: would overwrite input file\n" );
      if( option->log )
      {
        writeLog( "SCF file not written: would overwrite input file\n" );
      }
    }
  } /* if option->writeScf */


  /*
  ** Write phd file, if requested.
  */
  if( option->writePhd )
  {
    /*
    ** Construct sequence filename.
    */
    if( option->writePhdName == NULL )
    {
      sprintf( filename,
               "%.*s.phd.1",
               PHRED_PATH_MAX,
               fileRootName );
    }
    else
    {
      sprintf( filename,
               "%.*s",
               PHRED_PATH_MAX,
               option->writePhdName );
    }

    /*
    ** Add directory prefix, if requested.
    */
    if( option->pdOption == 1 )
    {
      /*
      ** Linux fails to write filename into filename with sprintf
      */
      pstrcpy( string, filename, PHRED_PATH_MAX );
      sprintf( filename,
               "%.*s%c%.*s",
               PHRED_PATH_MAX,
               option->phdDirName,
               PATHSEP,
               PHRED_PATH_MAX,
               string );
    }

    /*
    ** write phd file
    */
    if( writePhd( filename, seqName, option, phredData ) == ERROR )
    {
      return( ERROR );
    }
  } /* if option->writePhd */

	//
	// MPS: we should write our data here
	//

	





  /*
  ** Write polymorphism file, if requested.
  */
  if( option->writePolyData )
  {
    /*
    ** Construct polymorphism filename.
    */
    if( option->writePolyDataName == NULL )
    {
      sprintf( filename,
               "%.*s.poly",
               PHRED_PATH_MAX,
               fileRootName );
    }
    else
    {
      sprintf( filename,
               "%.*s",
               PHRED_PATH_MAX,
               option->writePolyDataName );
    }

    /*
    ** Add directory prefix, if requested.
    */
    if( option->ddOption == 1 )
	{
      /*
      ** Linux fails to write filename into filename with sprintf
      */
      pstrcpy( string, filename, PHRED_MAX_STRING_LEN );
      sprintf( filename,
               "%.*s%c%.*s",
               PHRED_PATH_MAX,
               option->polyDataDirName,
               PATHSEP,
               PHRED_PATH_MAX,
               string );
    }
    if( writePolyData( filename, seqName,
                       phredData->numPoint,
                       &(phredData->polyData) ) == ERROR )
    {
      return( ERROR );
	}

  }

  return( OK );
}


