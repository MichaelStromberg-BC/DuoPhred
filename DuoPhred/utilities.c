/** utilities.c **/

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
**    * utilities.c                                                    *     **
**    * benefits from ideas in code written by LaDeana Hillier and     *     **
**    * Tim Gleeson.                                                   *     **
**                                                                           **
*******************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "phred.h"

#ifdef ANSI_C
char *ourMalloc( int numByte )
#else
char *ourMalloc( numByte )
int numByte;
#endif
{
  char *cptr;

  if( numByte == 0 )
  {
    return( NULL );
  }

  if( ( cptr = (char *)malloc( numByte ) ) == NULL )
  {
    fprintf( stderr, "ourMalloc: error: unable to allocate memory\n" );
    return( NULL );
  }
  return( cptr );
}

#ifdef ANSI_C
int ourFree( char *cptr )
#else
int ourFree( cptr )
char *cptr;
#endif
{
  free( cptr );
  return( OK );
}


#ifdef ANSI_C
int stringMatch( char *seq1, int n1, char *seq2, int n2, int nmiss, int *indices )
#else
int stringMatch( seq1, n1, seq2, n2, nmiss, indices )
char *seq1, *seq2;
int n1, n2, nmiss;
int *indices;
#endif
{
    int i, j, d;
    int istart, iend;
    int i_miss, n_match;
    int mtable[100][5];


    if( n1 - n2 > nmiss )
    {
      return( 0 );
    }

    /*
    ** d = j - i is the "offset" between the two sequences
    */
    n_match = 0;
    for( d = -nmiss; d <= n2 + nmiss - n1; d++ )
    {
        if( d < 0 )
        {
          i_miss = -d;
          istart = -d;
        }
        else
        {
          i_miss = 0;
          istart = 0;
        }
        if( d > n2 - n1 )
        {
            iend = n2 - d;
            i_miss += n1 + d - n2;
        }
        else
        {
          iend = n1;
        }
        for( i = istart, j = d + i; i < iend; i++, j++ )
        {
            if( seq1[i] != seq2[j] && ++i_miss > nmiss )
            {
              goto nextd;
            }
        }

        /* indices (assuming they start at
        ** 0) of starting nucleotide in the
        ** searched sequence
        */
        mtable[n_match][0] = d + istart;

        /* nucleotide position in the query
        ** sequence where match starts (assuming
        ** query index starts with 0)
        */
        mtable[n_match][1] = istart;

        /* number of nucleotides in the
        ** match
        */
        mtable[n_match][2] = iend - istart;

        if( mtable[n_match][2] == n1 )
        {
          indices[n_match] = mtable[n_match][0];
        }

        /*
        ** Number of mismatches.
        */
        mtable[n_match][3] = i_miss;

        /*
        ** Number of matches.
        */
        n_match++;
        if( n_match >= 100 )
        {
          return( n_match );
        }

    nextd:;
    }

    return( n_match );
}


#ifdef ANSI_C
char *getFileName( char *fullPathName )
#else
char *getFileName( fullPathName )
char *fullPathName;
#endif
{
  char *cptr;

  if( ( cptr = strrchr( fullPathName, PATHSEP ) ) != NULL )
  {
    return( cptr + 1 );
  }
  return( fullPathName );
}


#ifdef ANSI_C
char *getTime( void )
#else
char *getTime()
#endif
{
  int i;
  time_t now_t;
  static char timstr[256];
  now_t = time( NULL );
  strcpy( timstr, ctime( &now_t ) );
  for( i = 0; i < (int)strlen( timstr ); ++i )
  {
    if( timstr[i] == '\n' )
    {
      timstr[i] = '\0';
      return( timstr );
    }
  }
  return( timstr );
}


