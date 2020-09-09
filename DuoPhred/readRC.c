/** readRC.c **/

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
#include <ctype.h>
#include "phred.h"

/*
** read .phredrc file
*/

/*
** valid entries:
**
** trim: <enzyme seq>
** view: true | false
** call: true | false
** logFile:  true | false
** bottom: true | false
** mag: n
** seq:
** seq: <filename>
** catSeq: <filename>
** seqDir: <directory name>
** seqType: fasta | xbap
** qual:
** qual: <filename>
** catQual: <filename>
** qualDir: <directory name>
** chromat:
** chromat: <filename>
** chromatDir: <directory name>
**
*/


#ifdef ANSI_C
int toLowerCase( char *cptr );
#else
int toLowerCase();
#endif

#ifdef ANSI_C
int readRC( Option *option )
#else
int readRC( option )
Option *option;
#endif
{
  int lineNum;
  int status;
  char line[PHRED_PATH_MAX];
  char savLine[PHRED_PATH_MAX];
  static char enzName[PHRED_PATH_MAX];
  static char seqName[PHRED_PATH_MAX];
  static char seqDirName[PHRED_PATH_MAX];
  static char qualName[PHRED_PATH_MAX];
  static char qualDirName[PHRED_PATH_MAX];
  static char scfName[PHRED_PATH_MAX];
  static char scfDirName[PHRED_PATH_MAX];
  char *cptr;
  FILE *fp;

  status = OK;
  lineNum = 0;
  if( ( fp = fopen( ".phredrc", "r" ) ) == NULL )
  {
    return( status );
  }

  while( fgets( line, PHRED_PATH_MAX - 1, fp ) != NULL )
  {
    ++lineNum;
    pstrcpy( savLine, line, PHRED_PATH_MAX );
    if( ( cptr = strtok( line, " \t\n\r:" ) ) != NULL )
    {
      /*
      ** view
      */
      if( strcmp( cptr, "view" ) == 0 )
      {
        if( ( cptr = strtok( NULL, " \t\n\r:" ) ) == NULL )
        {
          fprintf( stderr,
                   "error in line %d of .phredrc: missing argument\n",
                   lineNum );
          status = ERROR;
          option->errorFlag = 1;
          continue;
        }
        toLowerCase( cptr );
        if( strcmp( cptr, "false" ) == 0 )
        {
          option->edit = 0;
        }
        else
        if( strcmp( cptr, "true" ) == 0 )
        {
          option->edit = 1;
        }
        else
        {
          fprintf( stderr,
                   "error in line %d of .phredrc: bad value %s\n",
                   lineNum, cptr );
          status = ERROR;
          option->errorFlag = 1;
          continue;
        }
      }
      else
      /*
      ** trim
      */
      if( strcmp( cptr, "trim" ) == 0 )
      {
        if( ( cptr = strtok( NULL, " \t\n\r:" ) ) == NULL )
        {
          fprintf( stderr,
                   "error in line %d of .phredrc: missing argument\n",
                   lineNum );
          status = ERROR;
          option->errorFlag = 1;
          continue;
        }
        option->trim = 1;
        pstrcpy( enzName, cptr, PHRED_PATH_MAX );
        option->enzName = enzName;
      }
      else
      /*
      ** call
      */
      if( strcmp( cptr, "call" ) == 0 )
      {
        if( ( cptr = strtok( NULL, " \t\n\r:" ) ) == NULL )
        {
          fprintf( stderr,
                   "error in line %d of .phredrc: missing argument\n",
                   lineNum );
          status = ERROR;
          option->errorFlag = 1;
          continue;
        }
        toLowerCase( cptr );
        if( strcmp( cptr, "false" ) == 0 )
        {
          option->call = 0;
        }
        else
        if( strcmp( cptr, "true" ) == 0 )
        {
          option->call = 1;
        }
        else
        {
          fprintf( stderr,
                   "error in line %d of .phredrc: bad argument %s\n",
                   lineNum, cptr );
          status = ERROR;
          option->errorFlag = 1;
          continue;
        }
      }
      else
      /*
      ** seq
      */
      if( strcmp( cptr, "seq" ) == 0 )
      {
        if( ( cptr = strtok( NULL, " \t\n\r:" ) ) == NULL )
        {
          option->writeSeq = 1;
          option->sOption = 1;
          option->writeSeqName = NULL;
        }
        else
        {
          option->writeSeq = 1;
          option->sOption = 2;
          pstrcpy( seqName, cptr, PHRED_PATH_MAX );
          option->writeSeqName = seqName;
        }
      }
      else
      /*
      ** catSeq
      */
      if( strcmp( cptr, "catSeq" ) == 0 )
      {
        if( ( cptr = strtok( NULL, " \t\n\r:" ) ) == NULL )
        {
          fprintf( stderr,
                   "error in line %d of .phredrc: missing argument\n",
                   lineNum );
          status = ERROR;
          option->errorFlag = 1;
          continue;
        }
        option->writeSeq = 1;
        option->saOption = 1;
        pstrcpy( seqName, cptr, PHRED_PATH_MAX );
        option->writeSeqName = seqName;
      }
      else
      /*
      ** seqDir
      */
      if( strcmp( cptr, "seqDir" ) == 0 )
      {
        if( ( cptr = strtok( NULL, " \t\n\r:" ) ) == NULL )
        {
          fprintf( stderr,
                   "error in line %d of .phredrc: missing argument\n",
                   lineNum );
          status = ERROR;
          option->errorFlag = 1;
          continue;
        }
        option->writeSeq = 1;
        option->sdOption = 1;
        pstrcpy( seqDirName, cptr, PHRED_PATH_MAX );
        option->seqDirName = seqDirName;
      }
      else
      /*
      ** qual
      */
      if( strcmp( cptr, "qual" ) == 0 )
      {
        if( ( cptr = strtok( NULL, " \t\n\r:" ) ) == NULL )
        {
          option->writeQual = 1;
          option->qOption = 1;
          option->writeQualName = NULL;
        }
        else
        {
          option->writeQual = 1;
          option->qOption = 2;
          pstrcpy( qualName, cptr, PHRED_PATH_MAX );
          option->writeQualName = qualName;
        }
      }
      else
      /*
      ** catQual
      */
      if( strcmp( cptr, "catQual" ) == 0 )
      {
        if( ( cptr = strtok( NULL, " \t\n\r:" ) ) == NULL )
        {
          fprintf( stderr,
                   "error in line %d of .phredrc: missing argument\n",
                   lineNum );
          status = ERROR;
          option->errorFlag = 1;
          continue;
        }
        option->writeQual = 1;
        option->qaOption = 1;
        pstrcpy( qualName, cptr, PHRED_PATH_MAX );
        option->writeQualName = qualName;
      }
      else
      /*
      ** qualDir
      */
      if( strcmp( cptr, "qualDir" ) == 0 )
      {
        if( ( cptr = strtok( NULL, " \t\n\r:" ) ) == NULL )
        {
          fprintf( stderr,
                   "error in line %d of .phredrc: missing argument\n",
                   lineNum );
          status = ERROR;
          option->errorFlag = 1;
          continue;
        }
        option->writeQual = 1;
        option->qdOption = 1;
        pstrcpy( qualDirName, cptr, PHRED_PATH_MAX );
        option->qualDirName = qualDirName;
      }
      else
      /*
      ** chromat
      */
      if( strcmp( cptr, "chromat" ) == 0 )
      {
        option->writeScf = 1;
        if( ( cptr = strtok( NULL, " \t\n\r:" ) ) == NULL )
        {
          option->writeScfName = NULL;
        }
        else
        {
          pstrcpy( scfName, cptr, PHRED_PATH_MAX );
          option->writeScfName = scfName;
        }
      }
      else
      /*
      ** chromatDir
      */
      if( strcmp( cptr, "chromatDir" ) == 0 )
      {
        if( ( cptr = strtok( NULL, " \t\n\r:" ) ) == NULL )
        {
          fprintf( stderr,
                   "error in line %d of .phredrc: missing argument\n",
                   lineNum );
          status = ERROR;
          option->errorFlag = 1;
          continue;
        }
        option->writeScf = 1;
        option->scfDir = 1;
        pstrcpy( scfDirName, cptr, PHRED_PATH_MAX );
        option->scfDirName = scfDirName;
      }
      else
      /*
      ** logFile
      */
      if( strcmp( cptr, "logFile" ) == 0 )
      {
        if( ( cptr = strtok( NULL, " \t\n\r:" ) ) == NULL )
        {
          fprintf( stderr,
                   "error in line %d of .phredrc: missing argument\n",
                   lineNum );
          status = ERROR;
          option->errorFlag = 1;
          continue;
        }
        toLowerCase( cptr );
        if( strcmp( cptr, "false" ) == 0 )
        {
          option->log = 0;
        }
        else
        if( strcmp( cptr, "true" ) == 0 )
        {
          option->log = 1;
        }
        else
        {
          fprintf( stderr,
                   "error in line %d of .phredrc: bad argument %s\n",
                   lineNum, cptr );
          status = ERROR;
          option->errorFlag = 1;
          continue;
        }
      }
      else
      /*
      ** diag
      */
      /*
      ** mag
      */
      if( strcmp( cptr, "mag" ) == 0 )
      {
        if( ( cptr = strtok( NULL, " \t\n\r:" ) ) == NULL )
        {
          fprintf( stderr,
                   "error in line %d of .phredrc: missing argument\n",
                   lineNum );
          status = ERROR;
          option->errorFlag = 1;
          continue;
        }
        option->magnification = atoi( cptr );

      }
      else
      /*
      ** bottom
      */
      if( strcmp( cptr, "bottom" ) == 0 )
      {
        if( ( cptr = strtok( NULL, " \t\n\r:" ) ) == NULL )
        {
          fprintf( stderr,
                   "error in line %d of .phredrc: missing argument\n",
                   lineNum );
          status = ERROR;
          option->errorFlag = 1;
          continue;
        }
        toLowerCase( cptr );
        if( strcmp( cptr, "false" ) == 0 )
        {
          option->bottom = 0;
        }
        else
        if( strcmp( cptr, "true" ) == 0 )
        {
          option->bottom = 1;
        }
        else
        {
          fprintf( stderr,
                   "error in line %d of .phredrc: bad argument %s\n",
                   lineNum, cptr );
          status = ERROR;
          option->errorFlag = 1;
          continue;
        }
      }
      else
      /*
      ** seqType
      */
      if( strcmp( cptr, "seqType" ) == 0 )
      {
        if( ( cptr = strtok( NULL, " \t\n\r:" ) ) == NULL )
        {
          fprintf( stderr,
                   "error in line %d of .phredrc: missing argument\n",
                   lineNum );
          status = ERROR;
          option->errorFlag = 1;
          continue;
        }
        toLowerCase( cptr );
        if( strcmp( cptr, "fasta" ) == 0 )
        {
          option->seqType = 0;
        }
        else
        if( strcmp( cptr, "xbap" ) == 0 )
        {
          option->seqType = 1;
        }
        else
        {
          fprintf( stderr,
                   "error in line %d of .phredrc: bad argument %s\n",
                   lineNum, cptr );
          status = ERROR;
          option->errorFlag = 1;
          continue;
        }
      }
/*
** add entries immediately above
*/
      else
      if( cptr[0] == '#' )
      {
        continue;
      }
      else
      {
        fprintf( stderr,
                 "error in line %d of .phredrc: not an option %s\n",
                 lineNum, cptr );
        status = ERROR;
        option->errorFlag = 1;
        continue;
      }
    } /* if cptr */
  } /* while fgets */

  /*
  ** Close file.
  */
  fclose( fp );
  
  return( status );
}


#ifdef ANSI_C
int toLowerCase( char *cptr )
#else
int toLowerCase( cptr )
char *cptr;
#endif
{
  int i;
  int len;

  len = strlen( cptr );
  for( i = 0; i < len; ++i )
  {
    cptr[i] = tolower( cptr[i] );
  }

  return( OK );
}
