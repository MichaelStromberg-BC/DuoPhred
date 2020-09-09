/** autoPhred.c **/

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
#include <limits.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "phred.h"

#ifdef ANSI_C
int testCompress( char *ifnm, char *ofnm );
int uncompressFile( char *ifnm, char *ofnm, Option *option, int type );
int rmFile( char *filename );
#else
int testCompress();
int uncompressFile();
int rmFile();
#endif

#ifdef ANSI_C
int autoPhred( Option *option )
#else
int autoPhred( option )
Option *option;
#endif
{
  int ifil, nfil;
  int inType;
  int nProc;
  int numPoint;
  int numBase;
  int compressType;
  int status;
  int ret_status;
  char message[2*PHRED_PATH_MAX];
  char prcFileName[PHRED_PATH_MAX];
  char tmpFileName[PHRED_PATH_MAX];
  char *inFileName;
  ChromatData *chromatData;
  PhredData *phredData;
  struct STAT_STRUCT statBuffer;

  ret_status = OK;

  nfil = option->numInFile;
  inType = option->inType;

  /*
  ** Initialize histogram.
  */
  if( option->qualReport == 1 )
  {
    initReport();
  }

  if( option->log )
  {
    writeLog( "Processing:\n" );
  }

  /*
  ** Loop through input files.
  */
  nProc = 0;
  for( ifil = 0; ifil < nfil; ++ifil )
  {

    option->curInFile = ifil;
    inFileName = option->inFileName[ifil];

    /*
    ** Uncompress file if necessary.
    */
    if( ( compressType = testCompress( inFileName, tmpFileName ) ) !=
        NO_COMPRESS )
    {
      if( uncompressFile( tmpFileName, prcFileName, option, compressType ) ==
          ERROR )
      {
        pstrcpy( prcFileName, inFileName, PHRED_PATH_MAX );
      }
    }
    else
    {
      pstrcpy( prcFileName, inFileName, PHRED_PATH_MAX );
    }

    /*
    ** Read data file.
    */
    chromatData = readData( prcFileName, &status );

    /*
    ** Delete temporary uncompressed file, if necessary.
    */
    if( compressType != NO_COMPRESS )
    {
      rmFile( prcFileName );
    }

    /*
    ** Print file reading status to standard out and standard error.
    ** Exit program on fatal error.
    */
    if( reportFileStatus( option, inFileName, status ) == ERROR )
    {
      if( chromatData != NULL )
      {
        freeChromatData( chromatData );
      }
      return( ERROR );
    }

    /*
    ** Skip if chromatData is NULL: indicates unable to read file.
    */
    if( chromatData == NULL )
    {
      continue;
    }

    /*
    ** Sanity test.
    */
    if( ( chromatData->numPoint      == 0 ||
          chromatData->maxTraceValue == 0.0 ) &&
        status == 0 )
    {
      if( option->tagOption )
      {
/*
** OK
*/
        sprintf( message,
                 "FATAL_ERROR: %.*s: inconsistent internal status: numPoint: %d  maxTaceVale: %f  status: %d\n",
                 PHRED_PATH_MAX,
                 inFileName,
                 numPoint,
                 chromatData->maxTraceValue,
                 status );
      }
      else
      {
/*
** OK
*/
        sprintf( message,
                 "    fatal error while processing file %.*s: inconsistent internal status: numPoint: %d  maxTaceVale: %f  status: %d\n",
                 PHRED_PATH_MAX,
                 inFileName,
                 numPoint,
                 chromatData->maxTraceValue,
                 status );
      }
      fprintf( stderr, "%s", message );

      /*
      ** log error if requested
      */
      if( option->log )
      {
        writeLog( message );
      }
      return( ERROR );
    }

    /*
    ** Get input file inode.
    */
    FILE_STATUS( inFileName, &statBuffer  );
    option->inInode = statBuffer.st_ino;

    option->readFileName = prcFileName;

    /*
    ** Create phredData structure and
    ** copy necessary data to it.
    */
    phredData = chromat2phred( chromatData, option );

    /*
    ** Free chromat data.
    */
    freeChromatData( chromatData );

    numPoint = phredData->numPoint;
    numBase  = phredData->numBase[FIL];

    /*
    ** Count number of files processed.
    */
    if( numPoint > 0 &&
        phredData->maxTraceValue != 0.0 )
    {
      ++nProc;
    }

    /*
    ** Set data set 'pointers'.
    */
    if( option->call == 1 )
    {
      phredData->dataSet[0] = LCL;
      phredData->dataSet[1] = FIL;
    }
    else
    {
      phredData->dataSet[0] = FIL;
      phredData->dataSet[1] = LCL;
    }

    /*
    ** Original trimming requires chromat base calls to locate peaks.
    */
    if( option->trim == 1 )
    {
      if( numPoint > 0 &&
          phredData->numBase[FIL] > 0 &&
          phredData->maxTraceValue != 0.0 )
      {
        /*
        ** Autoclip.
        */
        if( trimSeq( phredData, option ) == ERROR )
        {
          if( option->tagOption )
          {
/*
** OK
*/
            sprintf( message,
                     "PROCESSING_ERROR: %.*s: trimming error\n",
                     PHRED_PATH_MAX,
                     inFileName );
          }
          else
          {
/*
** OK
*/
            sprintf( message,
                     "    trimming error in file %.*s\n",
                     PHRED_PATH_MAX,
                     inFileName );
          }
          fprintf( stderr, "%s\n", message );
          if( option->log )
          {
            writeLog( message );
          }
        }
      }
      else
      {
        if( numPoint > 0 && phredData->numBase[FIL] == 0 )
        {
          if( option->tagOption )
          {
/*
** OK
*/
            sprintf( message,
                     "PROCESSING_ERROR: %.*s:  trimming error: '-trim' requires chromat trace base calls\n",
                     PHRED_PATH_MAX,
                     inFileName );
          }
          else
          {
/*
** OK
*/
            sprintf( message,
                     "    trimming error in file %.*s: '-trim' requires chromat base calls\n",
                     PHRED_PATH_MAX,
                     inFileName );
          }
          fprintf( stderr, "%s", message );
          if( option->log )
          {
            writeLog( message );
          }
        }
      }
    }

    /*
    ** Call bases if requested.
    */
    if( option->call )
    {
      if( numPoint > 0 &&
          phredData->maxTraceValue != 0.0 )
      {
        /*
        ** Call bases if trace exists. User was warned earlier
        ** of missing trace.
        */
        if( callSeq( inFileName, phredData, option, &status ) != ERROR )
        {
          /*
          ** Set trim points of local sequence.
          */
          if( option->trim == 1 && phredData->numBase[FIL] > 0 ) 
          {
            trimSet( phredData );
          }
          else
          if( option->trim == 2 )
          {
            trimAlt( phredData, option->enzName );
          }

          /*
          ** Trimming information for PHD file.
          */
          trimPhred( phredData );
        }
      }
    }
    else
    {
      /*
      ** It's possible to trim sequence from the bases in the chromat file.
      */
      if( numPoint > 0 &&
          phredData->maxTraceValue != 0.0 &&
          numBase > 0 )
      {
        if( option->trim == 1 )
        {
          trimSet( phredData );
        }
      }
    }

    /*
    ** Evaluate trace quality for read histogram.
    */
    if( option->qualReport == 1 &&
        numPoint > 0 &&
        phredData->maxTraceValue != 0.0 &&
        makeReport( phredData ) == ERROR )
    {
      if( option->tagOption )
      {
/*
** OK
*/
        sprintf( message,
                 "PROCESSING_ERROR: %.*s: error making histogram\n",
                 PHRED_PATH_MAX,
                 inFileName );
      }
      else
      {
/*
** OK
*/
        sprintf( message,
                 "    error making histogram %.*s\n",
                 PHRED_PATH_MAX,
                 inFileName );
      }
      fprintf( stderr, "%s\n", message );
      if( option->log )
      {
        writeLog( message );
      }
    }

    /*
    ** Write data.
    */
    if( writeData( phredData, prcFileName, option ) == ERROR )
    {
      if( option->tagOption )
      {
/*
** OK
*/
        sprintf( message,
                 "FILE_ERROR: %.*s: write error\n",
                 PHRED_PATH_MAX,
                 prcFileName );
      }
      else
      {
/*
** OK
*/
        sprintf( message,
                 "    write error in file %.*s\n",
                 PHRED_PATH_MAX,
                 prcFileName );
      }
      fprintf( stderr, "%s\n", message );
      if( option->log )
      {
        writeLog( message );
      }

      if( option->writePolyData == 1 )
      {
        freePolyData( &(phredData->polyData) );
      }

      freePhredData( phredData );

      ret_status = ERROR;

      continue;
    }

    if( option->writePolyData == 1 )
    {
      freePolyData( &(phredData->polyData) );
    }

    /*
    ** Free phredData memory.
    */
    freePhredData( phredData );

  }

  if( option->log )
  {
    writeLog( "\n\nCompleted processing.\n" );
  }

  /*
  ** Write histogram.
  */
  if( option->qualReport == 1 && nProc > 0 )
  {
    if( writeReport( option ) == ERROR )
    {
      if( option->tagOption )
      {
/*
** OK
*/
        sprintf( message,
                 "FILE_ERROR: %.*s: report write error\n",
                 PHRED_PATH_MAX,
                 option->qualReportFileName );
      }
      else
      {
/*
** OK
*/
        sprintf( message, "    report write error\n" );
      }
      fprintf( stderr, "%s\n", message );
      if( option->log )
      {
        writeLog( message );
      }
    }
  }

  return( ret_status );
}

