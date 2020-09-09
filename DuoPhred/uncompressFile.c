/** uncompressFile.c **/

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
int uncompressFile( char *ifnm, char *ofnm, Option *option, int type )
#else
int uncompressFile( ifnm, ofnm, option, type )
char *ifnm;
char *ofnm;
Option *option;
int type;
#endif
{
  int i, len;
  char exeName[2*PHRED_PATH_MAX];
  char tmpName[2*PHRED_PATH_MAX];
  char string[4*PHRED_PATH_MAX];
  char *cptr;

  /*
  ** Build (path) name of compression executable.
  */
  if( option->filCompressExeDirOption )
  {
    if( type == GZ_COMPRESS )
    {
      sprintf( exeName, "%.*s%c%s",
               PHRED_PATH_MAX,
               option->filCompressExeDir,
               PATHSEP,
               "gunzip" );
    }
    else
    if( type == Z_COMPRESS )
    {
      sprintf( exeName, "%.*s%c%s",
               PHRED_PATH_MAX,
               option->filCompressExeDir,
               PATHSEP,
               "uncompress" );
    }
    else
    if( type == BZ_COMPRESS )
    {
      sprintf( exeName, "%.*s%c%s",
               PHRED_PATH_MAX,
               option->filCompressExeDir,
               PATHSEP,
               "bzip2 -d -k" );
    }
    else
    {
      sprintf( exeName, "%.*s%c%s", 
               PHRED_PATH_MAX,
               option->filCompressExeDir,
               PATHSEP,
               "uncompress" );
    }
  }
  else
  {
    if( type == GZ_COMPRESS )
    {
      strcpy( exeName, "gunzip" );
    }
    else
    if( type == Z_COMPRESS )
    {
      strcpy( exeName, "uncompress" );
    }
    else
    if( type == BZ_COMPRESS )
    {
      strcpy( exeName, "bzip2 -d -k" );
    }
    else
    {
      strcpy( exeName, "uncompress" );
    }
  }

  /*
  ** Build path of temporary directory.
  */
  if( option->filCompressTmpDirOption )
  {
    sprintf( tmpName,
             "%.*s",
             PHRED_PATH_MAX,
             option->filCompressTmpDir );
  }
  else
  {
    sprintf( tmpName,
             "%s",
             "/usr/tmp" );
  }

  /*
  ** Build name of uncompressed file.
  */
  if( ( cptr = strrchr( ifnm, '/' ) ) == NULL )
  {
    cptr = ifnm;
  }
  else
  {
    ++cptr;
  }
  sprintf( string,
           "%.*s/%.*s",
           PHRED_PATH_MAX - 8,
           tmpName,
           PHRED_PATH_MAX - 8,
           cptr );
  len = strlen( string );
  if( type == GZ_COMPRESS &&
      string[len-1] == 'z' &&
      string[len-2] == 'g' &&
      string[len-3] == '.' )
  {
    len -= 3;
  }
  else
  if( type == Z_COMPRESS &&
      string[len-1] =='Z' &&
      string[len-2] == '.' )
  {
    len -= 2;
  }
  else
  if( type == BZ_COMPRESS &&
      string[len-1] == '2' &&
      string[len-2] == 'z' &&
      string[len-3] == 'b' &&
      string[len-4] == '.' )
  {
    len -= 4;
  }
/*
  else
  {
    len -= 2;
  }
*/

  for( i = 0; i < len; ++i )
  {
    ofnm[i] = string[i];
  }
  ofnm[len] = '\0';

  /*
  ** Uncompress file in temporary directory.
  */
  sprintf( string,
           "%.*s -f -c %.*s > %.*s",
           PHRED_PATH_MAX,
           exeName,
           PHRED_PATH_MAX,
           ifnm,
           PHRED_PATH_MAX,
           ofnm );
  if( system( string ) != 0 )
  {
    fprintf( stderr,
             "uncompressFile: unable to uncompress file %s\n",
             ifnm );
    system( string );
    return( ERROR );
  }
  
  return( OK );
}
