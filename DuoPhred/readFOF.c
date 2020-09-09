/** readFOF.c **/

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

static char *sfnm = NULL;

#ifdef ANSI_C
int readFOF( int argc, char *filename, Option *option )
#else
int readFOF( argc, filename, option )
int argc;
char *filename;
Option *option;
#endif
{
  int i, j, n;
  char line[PHRED_PATH_MAX];
  char *cptr;
  FILE *fp;

  if( ( fp = fopen( filename, "r" ) ) == NULL )
  {
    fprintf( stderr, "unable to open file %s\n", filename );
    return( ERROR );
  }
  n = 0;
  while( fgets( line, PHRED_PATH_MAX, fp ) != 0 )
  {
    ++n;
  }
  if( ( option->inFileName = (char **)realloc( option->inFileName,
                     ( argc + n ) * sizeof( char * ) ) ) == NULL )
  {
    fprintf( stderr, "unable to allocate memory\n" );
    fclose( fp );
    return( ERROR );
  }
  if( ( sfnm = (char *)malloc( ( n + 1 ) * ( PHRED_PATH_MAX + 1 ) * sizeof( char ) ) ) == NULL )
  {
    fprintf( stderr, "unable to allocate memory\n" );
    fclose( fp );
    return( ERROR );
  }
  rewind( fp );
  cptr = sfnm;
  for( i = 0; i < n; ++i )
  {
    fgets( line, PHRED_PATH_MAX, fp );
    option->inFileName[option->numInFile] = cptr;
    for( j = 0; j < (int)strlen( line ); ++j )
    {
      if( line[j] == ' ' )
      {
        continue;
      }
      if( line[j] == '\n' ||
          line[j] == '\0' )
      {
        *cptr = '\0';
        /*
        ** retain pointer if some characters were saved
        */
        if( cptr != option->inFileName[option->numInFile] )
        {
          ++option->numInFile;
        }
        ++cptr;
        break;
      }
      *cptr = line[j];
      ++cptr;
    }
  }
  fclose( fp );

  return( OK );
}

#ifdef ANSI_C
int initFOF( void )
#else
int initFOF( )
#endif
{
  sfnm = NULL;
  return( OK );
}

#ifdef ANSI_C
int freeFOF( void )
#else
int freeFOF( )
#endif
{
  if( sfnm != NULL )
  {
    free( sfnm );
  }
  return( OK );
}
