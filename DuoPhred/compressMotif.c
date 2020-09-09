/** compressMotif.c **/

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

#define STRONG_3	3
#define STRONG_4	4

#define WEAK_3		13
#define WEAK_4		14
#define WEAK_5		15
#define WEAK_6		16

#ifdef ANSI_C
static int check_stem( int loop_strength, int stem_size, char *stem_seq )
#else
static int check_stem( loop_strength, stem_size, stem_seq )
int loop_strength;
int stem_size;
char *stem_seq;
#endif
{
  int i;
  int nGC;

  if( stem_size == 2 )
  {
    if( loop_strength <= STRONG_4 )
    {
      return( 1 );
    }
    else
    {
      return( 0 );
    }
  }
  else
  if( stem_size == 3 )
  {
    if( stem_seq[stem_size-3] == 'G' &&
        stem_seq[stem_size-2] == 'C' &&
        stem_seq[stem_size-1] == 'C' )
    {
      if( loop_strength <= WEAK_5 )
      {
        return( 1 );
      }
      else
      {
        return( 0 );
      }
    }
    else
    if( stem_seq[stem_size-3] == 'G' &&
        stem_seq[stem_size-2] == 'G' &&
        stem_seq[stem_size-1] == 'G' )
    {
      if( loop_strength <= WEAK_3 )
      {
        return( 1 );
      }
      else
      {
        return( 0 );
      }
    }
    else
    {
      if( loop_strength <= STRONG_4 )
      {
        return( 1 );
      }
      else
      {
        return( 0 );
      }
    }
  }
  else
  if( stem_size == 4 )
  {
    if( stem_seq[stem_size-3] == 'G' &&
        stem_seq[stem_size-2] == 'C' &&
        stem_seq[stem_size-1] == 'C' )
    {
      if( loop_strength == WEAK_4 )
      {
        return( 1 );
      }
    }

    nGC = 0;
    for( i = 0; i < stem_size && nGC < 3; ++i )
    {
      if( stem_seq[i] == 'G' ||
          stem_seq[i] == 'C' )
      {
        ++nGC;
      }
    }

    if( nGC >= 3 &&
        ( loop_strength <= WEAK_3 ) )
    {
      return( 1 );
    }
    else
    {
      return( 0 );
    }
  }
  else
  if( stem_size >= 5 )
  {
    nGC = 0; 
    for( i = 0; i < stem_size && nGC < 3; ++i )
    {
      if( stem_seq[i] == 'G' ||
          stem_seq[i] == 'C' )
      {
        ++nGC;
      }
    }

    if( nGC >= 3 &&
        ( loop_strength <= WEAK_6 ) )
    {
      return( 1 );
    }
    else
    {
      return( 0 );
    }
  }

  fprintf( stderr, "checkStem: error: unexpected stem sequence\n" );
    
  return( -1 );
}

#ifdef ANSI_C
int compressMotif( int numBase, char *seq, int *compress )
#else
int compressMotif( numBase, seq, compress )
int numBase;
char *seq;
int *compress;
#endif
{
  int i, j, k;
  int result;
  int numByte;
  int stem_size;
  int loop_size;
  int loop_strength;
  char comp[128];
  char *stem_seq;

  numByte = ( numBase + 1 ) * sizeof( char );
  stem_seq = (char *)malloc( numByte );
  if( stem_seq == NULL )
  {
    fprintf( stderr, "compressMotif: error: unable to allocate memory\n" );
    return( ERROR );
  }

  comp['A'] = 'T';
  comp['C'] = 'G';
  comp['G'] = 'C';
  comp['T'] = 'A';
  comp['N'] = 'N';
  comp['-'] = '-';

  for( i = 6; i < numBase; ++i )
  {
    for( loop_size = 3; loop_size <= 6; ++loop_size )
    {
      if( loop_size == 3 )
      {
        if( seq[i-3] == 'G' &&
            seq[i-1] == 'A' )
        {
          loop_strength = STRONG_3;
        }
        else
        {
          loop_strength = WEAK_3;
        }
      }
      else
      if( loop_size == 4 )
      {
        if( ( seq[i-4] == 'G' &&
              seq[i-1] == 'A' ) ||
            ( seq[i-4] == 'C' &&
              seq[i-1] == 'G' ) )
        {
          loop_strength = STRONG_4;
        }
        else
        {
          loop_strength = WEAK_4;
        }
      }
      else
      if( loop_size == 5 )
      {
        loop_strength = WEAK_5;
      }
      else
      if( loop_size == 6 )
      {
        loop_strength = WEAK_6;
      }
      else
      {
        fprintf( stderr, "compressMotif: error: undefined loop size\n" );
        free( stem_seq );
        return( ERROR );
      }

      k = i - loop_size - 1;
      stem_size = 0;

      for( j = i; j < numBase && k >= 0; ++j, --k )
      {
        stem_seq[stem_size] = comp[(int)seq[k]];
        ++stem_size;
        stem_seq[stem_size] = '\0';

        if( stem_size > 1 &&
            ( ( stem_seq[stem_size-2] == 'G' &&
                stem_seq[stem_size-1] == 'G' ) ||
              ( stem_seq[stem_size-2] == 'C' &&
                stem_seq[stem_size-1] == 'C' ) ) )
        {
          result = check_stem( loop_strength, stem_size, stem_seq );
          if( result == -1 )
          {
            free( stem_seq );
            return( ERROR );
          }
          if( result == 1 )
          {
            compress[j-1] = 1;
          }
        }

        if( seq[j] != comp[(int)seq[k]] )
        {
          break;
        }
      }
    }
  }

  /*
  ** Free memory.
  */
  free( stem_seq );

  return( OK );
}

