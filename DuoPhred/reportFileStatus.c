/** reportFileStatus.c **/

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
#include "phred.h"

#ifdef ANSI_C
int reportFileStatus( Option *option, char *inFileName, int statusCode )
#else
int reportFileStatus( option, inFileName, statusCode )
Option *option;
char *inFileName;
int statusCode;
#endif
{
  int iret;
  char message[2*PHRED_PATH_MAX];

  iret = OK;

  if( statusCode == 0 || statusCode == 3 )
  {
    if( option->tagOption )
    {
      sprintf( message,
               "PROCESS: %.*s\n",
               PHRED_PATH_MAX,
               inFileName );
    }
    else
    {
      sprintf( message,
               "  %.*s\n",
               PHRED_PATH_MAX,
               inFileName );
    }
    printf( "%s", message );

    /*
    ** log filename if requested
    */
    if( option->log )
    {
      writeLog( message );
    }
  }
  else
  if( statusCode == 1 )
  {
    if( option->tagOption )
    {
      sprintf( message,
               "FILE_ERROR: %.*s: file read error\n",
               PHRED_PATH_MAX,
               inFileName );
    }
    else
    {
      sprintf( message,
               "  read error in file %.*s\n",
               PHRED_PATH_MAX,
               inFileName );
    }
    fprintf( stderr, "%s", message );

    /*
    ** log error if requested
    */
    if( option->log )
    {
      writeLog( message );
    }
  }
  else
  if( statusCode == 2 )
  {
    if( option->tagOption )
    {
      sprintf( message,
               "FILE_ERROR: %.*s: trace data missing\n",
               PHRED_PATH_MAX,
               inFileName );
    }
    else
    {
      sprintf( message,
               "  %.*s: trace data missing\n",
               PHRED_PATH_MAX,
               inFileName );
    }
    fprintf( stderr, "%s", message );

    /*
    ** log error if requested
    */
    if( option->log )
    {
      writeLog( message );
    }
  }
  else
  if( statusCode == -1 )
  {
    if( option->tagOption )
    {
      sprintf( message,
               "FATAL_ERROR: %.*s: error while reading\n",
               PHRED_PATH_MAX,
               inFileName );
    }
    else
    {
      sprintf( message,
               "  fatal error while reading %.*s\n",
               PHRED_PATH_MAX,
               inFileName );
    }
    fprintf( stderr, "%s", message );

    /*
    ** log error if requested
    */
    if( option->log )
    {
      writeLog( message );
    }
    iret = ERROR;
  }
  else
  {
    if( option->tagOption )
    {
      sprintf( message,
               "FATAL_ERROR: %.*s: unknown status code\n",
               PHRED_PATH_MAX,
               inFileName );
    }
    else
    {
      sprintf( message,
               "  %.*s: fatal error: unknown status code while reading\n",
               PHRED_PATH_MAX,
               inFileName );
    }
    fprintf( stderr, "%s", message );

    /*
    ** log error if requested
    */
    if( option->log )
    {
      writeLog( message );
    }
    iret = ERROR;
  }

  return( iret );
}
