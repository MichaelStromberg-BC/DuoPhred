/** logFile.c **/

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
#include <time.h>
#ifndef _WIN32
#include <sys/time.h>
#else
#include <time.h>
#endif
#include "phred.h"

/*
** log file functions
*/
static FILE *sfp = NULL;


#ifdef ANSI_C
int writeLog( char *string );
#else
int writeLog();
#endif


/*
** log parameters
*/
#ifdef ANSI_C
int logParam( int argc, char **argv )
#else
int logParam( argc, argv )
int argc;
char **argv;
#endif
{
  int i;

  writeLog( "Command line: " );
  for( i = 1; i < argc; ++i )
  {
    writeLog( argv[i] );
    writeLog( " " );
  }
  writeLog( "\n\n" );

  return( OK );
}



/*
** open logfile and write time/date
*/
#ifdef ANSI_C
int openLog( void )
#else
int openLog()
#endif
{
  time_t pTime;

  if( ( sfp = fopen( "phred.log", "a" ) ) == NULL )
  {
    fprintf( stderr, "unable to open log file\n" );
    return( ERROR );
  }
  time( &pTime );
  fprintf( sfp, "\n\n---------- PHRED log : %s \n",
                ctime( &pTime ) );
  return( OK );
}


/*
** write string to log file
*/
#ifdef ANSI_C
int writeLog( char *string )
#else
int writeLog( string )
char *string;
#endif
{
  fprintf( sfp, "%s", string );
  return( OK );
}



/*
** close log file
*/
#ifdef ANSI_C
int closeLog( void )
#else
int closeLog()
#endif
{
  fclose( sfp );
  return( OK );
}
