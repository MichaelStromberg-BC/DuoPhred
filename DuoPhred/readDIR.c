/** readDIR.c **/

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
#include <sys/types.h>
#if defined(_WIN32)
#include <io.h>
#else
#include <dirent.h>
#endif
#include "phred.h"

static char *sdnm = NULL;


#if defined(_WIN32)

/*
** Windows 32 version
**
*/
#ifdef ANSI_C
int readDIR( int argc, char *dirname, Option *option )
#else
int readDIR( argc, dirname, option )
int argc;
char *dirname;
Option *option;
#endif
{
  int i, n;
  char line[3*PHRED_PATH_MAX];
  char *cptr;
  long hFile;
  static struct _finddata_t c_file;
  char pathpattern[2*PHRED_PATH_MAX];

  sprintf( pathpattern,
           "%.*s/%s",
           PHRED_PATH_MAX,
           dirname,
           "*.*" );

  if( (hFile = _findfirst( pathpattern, &c_file )) == -1L )
  {
    fprintf( stderr, "no files found\n");
    return( ERROR );
  }
  else
  {
    n = 1;
    while( _findnext( hFile, &c_file ) == 0 )
    n++;
  }

  _findclose( hFile );
  hFile = 0;

  if( ( option->inFileName = (char **)realloc( option->inFileName,
                             ( argc + n ) * sizeof( char * ) ) ) == NULL )
  {
    fprintf( stderr, "unable to allocate memory\n" );
    return( ERROR );
  }

  if( ( sdnm = (char *)malloc( ( n + 1 ) * ( PHRED_PATH_MAX + 1 ) * sizeof( char ) ) ) == NULL )
  {
    fprintf( stderr, "unable to allocate memory\n" );
    return( ERROR );
  }

  cptr = sdnm;
  for( i = 0; i < n; ++i )
  {
    if(hFile)
    {
      if(_findnext( hFile, &c_file ) != 0)
      continue;
    }
    else
    {
      hFile = _findfirst( pathpattern, &c_file );
    }

    if( strcmp( c_file.name, "."  ) == 0 ||
        strcmp( c_file.name, ".." ) == 0 )
    {
      continue;
    }

    sprintf( line,
             "%.*s%c%.*s",
             PHRED_PATH_MAX,
             dirname,
             PATHSEP,
             PHRED_PATH_MAX,
             c_file.name );
    strncpy( cptr, line, strlen( line ) + 1 );
    option->inFileName[option->numInFile] = cptr;
    ++option->numInFile;
    cptr += strlen( line ) + 1;
  }

  return( OK );
}




#else




/*
** UNIX version
**
*/
#ifdef ANSI_C
int readDIR( int argc, char *dirname, Option *option )
#else
int readDIR( argc, dirname, option )
int argc;
char *dirname;
Option *option;
#endif
{
  int i, n;
  char line[3*PHRED_PATH_MAX];
  char *cptr;
  DIR *dp;
  struct dirent *dep;

  if( ( dp = opendir( dirname ) ) == NULL )
  {
    fprintf( stderr, "unable to open directory\n" );
    return( ERROR );
  }

  n = 0;
  while( ( dep = readdir( dp ) ) != NULL )
  {
    /* dep->d_name */
    ++n;
  }

  if( ( option->inFileName = (char **)realloc( option->inFileName,
                     ( argc + n ) * sizeof( char * ) ) ) == NULL )
  {
    fprintf( stderr, "unable to allocate memory\n" );
    closedir( dp );
    return( ERROR );
  }
  if( ( sdnm = (char *)malloc( ( n + 1 ) * ( PHRED_PATH_MAX + 1 ) * sizeof( char ) ) ) == NULL )
  {
    fprintf( stderr, "unable to allocate memory\n" );
    closedir( dp );
    return( ERROR );
  }
  rewinddir( dp );
  cptr = sdnm;
  for( i = 0; i < n; ++i )
  {
    if( ( dep = readdir( dp ) ) == NULL )
    {
      continue;
    }
    if( strcmp( dep->d_name, "."  ) == 0 ||
        strcmp( dep->d_name, ".." ) == 0 )
    {
      continue;
    }
    sprintf( line,
             "%.*s%c%.*s",
            PHRED_PATH_MAX,
            dirname,
            PATHSEP,
            PHRED_PATH_MAX,
            dep->d_name );
    strncpy( cptr, line, strlen( line ) + 1 );
    option->inFileName[option->numInFile] = cptr;
    ++option->numInFile;
    cptr += strlen( line ) + 1;

  }

  closedir( dp );

  return( OK );
}




#endif


#ifdef ANSI_C
int initDIR( void )
#else
int initDIR( )
#endif
{
  sdnm = NULL;
  return( OK );
}



#ifdef ANSI_C
int freeDIR( void )
#else
int freeDIR( )
#endif
{
  if( sdnm != NULL )
  {
    free( sdnm );
  }
  return( OK );
}
