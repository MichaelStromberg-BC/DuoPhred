/** trimSeq.c **/

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
**    * trimSeq.c was written by LaDeana Hillier.                       *    **
**                                                                           **
*******************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "phred.h"

#define SIDEBAND_CUTOFF 0.8
#define NONCALLED_OVER_CALLED_CUTOFF 0.25
#define OVERALL_TRACE_QUAL_CUTOFF 0.27
#define STEP_SIZE 8
#define LAST_ALLOWED_BASE 550

#ifdef ANSI_C
FLOAT get_area( FLOAT *trace, int startp, int endp )
#else
FLOAT get_area( trace, startp, endp )
FLOAT *trace;
int startp;
int endp;
#endif
{
  int i;
  FLOAT sum=0.0;
  for( i = startp; i < endp; i++ )
  {
    sum += trace[i];
  }

  return( sum );
}

#ifdef ANSI_C
FLOAT max_area( FLOAT *tx, FLOAT *ty, FLOAT *tz, int stp, int endp )
#else
FLOAT max_area( tx, ty, tz, stp, endp )
FLOAT *tx;
FLOAT *ty;
FLOAT *tz;
int stp;
int endp;
#endif
{
  FLOAT x, y, z;

  x = get_area( tx, stp, endp );
  y = get_area( ty, stp, endp );
  z = get_area( tz, stp, endp );

  if( x > y )
  {
    if( z > x )
    {
      return( z );
    }
    else
    {
      return( x );
    }
  }
  else
  {
    if( z > y )
    {
      return( z );
    }
    else
    {
      return( y );
    }
  }
}

/*
** Returns the position half way between base and the following base.
*/
#ifdef ANSI_C
static FLOAT one_half_forwards( PhredData *phredData, int base )
#else
static FLOAT one_half_forwards( phredData, base)
PhredData *phredData;
int base;
#endif
{
  FLOAT pos;

  if( ( base + 1 ) < phredData->numBase[FIL] )
  {
    pos = ( phredData->baseLoc[FIL][base] + phredData->baseLoc[FIL][base+1] ) / 2.0;
  }
  else
  {
    /*
    ** Last base is a special case. We should guestimate.
    */
    pos = phredData->baseLoc[FIL][base] +
          ( phredData->baseLoc[FIL][base] - phredData->baseLoc[FIL][base-1] ) / 2.0;
    if( pos >= phredData->numPoint )
    {
      pos = (FLOAT)( phredData->numPoint - 1 );
    }
  }

  return( pos );
}


/*
** Returns the position half way between base and the preceding base.
*/
#ifdef ANSI_C
static FLOAT one_half_backwards( PhredData *phredData, int base )
#else
static FLOAT one_half_backwards( phredData, base )
PhredData *phredData;
int base;
#endif
{
  FLOAT pos;

  if( base > 0 )
  {
    pos = phredData->baseLoc[FIL][base] -
          ( ( phredData->baseLoc[FIL][base] - phredData->baseLoc[FIL][base-1] ) ) / 2;
  }
  else
  {
    /*
    ** Last base is a special case. We should guestimate.
    */
    pos = phredData->baseLoc[FIL][base] -
          ( phredData->baseLoc[FIL][base+1] - phredData->baseLoc[FIL][base] ) / 2.0;
    if( pos < 0.0 )
    {
      pos = 0.0;
    }
  }

  return( pos );

}


#ifdef ANSI_C
static void  SeqQual_nonCalledOverCalled( PhredData *phredData )
#else
static void  SeqQual_nonCalledOverCalled( phredData )
PhredData *phredData;
#endif
{
  int i, j;
  int pos, start_pos, end_pos;
  int half_window_size = 8;  /* look at three bases on either side */
  FLOAT area;
  int end_sum, start_sum;

  for( i = 0; i < phredData->numBase[FIL]; i++ )
  {
    phredData->qualIndex[i] = 0.0;
  }

  phredData->qualType = 15;

  for( j = 0; j < phredData->numBase[FIL]; j++ )
  {

    end_sum = j + half_window_size;
    start_sum = j - half_window_size;

    if( end_sum >= phredData->numBase[FIL] )
    {
      end_sum = ( phredData->numBase[FIL] ) - 1;
    }
    if( start_sum < 0 )
    {
      start_sum = 0;
    }

    for( i = start_sum; i <= end_sum; i++ )
    {
      pos = phredData->baseLoc[FIL][i];

      start_pos = one_half_backwards( phredData, i );
      end_pos = one_half_forwards( phredData, i );

      /*
      ** bge
      */
      if( end_pos <= start_pos )
      {
        continue;
      }

      switch( phredData->base[FIL][i] )
      {
        case 'A':
        case 'a':
          area = get_area( phredData->trace[0], start_pos, end_pos );
          if( area > 0.0 )
          {
            phredData->qualIndex[j] += max_area( phredData->trace[1],
                                                 phredData->trace[2],
                                                 phredData->trace[3],
                                                 start_pos, end_pos ) / area;
          }
          break;
        case 'C':
        case 'c':
          area = get_area( phredData->trace[1], start_pos, end_pos );
          if( area > 0.0 )
          {
            phredData->qualIndex[j] += max_area( phredData->trace[0],
                                                 phredData->trace[2],
                                                 phredData->trace[3],
                                                 start_pos, end_pos ) / area;
          }
          break;
        case 'G':
        case 'g':
          area = get_area( phredData->trace[2], start_pos, end_pos );
          if( area > 0.0 )
          {
            phredData->qualIndex[j] += max_area( phredData->trace[1],
                                                 phredData->trace[0],
                                                 phredData->trace[3],
                                                 start_pos, end_pos ) / area;
          }
          break;
        case 'T':
        case 't':
          area = get_area( phredData->trace[3], start_pos, end_pos );
          if( area > 0.0 )
          {
            phredData->qualIndex[j] += max_area( phredData->trace[1],
                                                 phredData->trace[2],
                                                 phredData->trace[0],
                                                 start_pos, end_pos ) / area;
          }
          break;
        default:
            phredData->qualIndex[j] += 1.0;
          break;
      }
    }
            phredData->qualIndex[j] = phredData->qualIndex[j] / ( end_sum - start_sum + 1 );
  }

   /*
   ** identify quality index as type 15
   */
   phredData->qualType = 15;
}

  /*
   * MODULE     SeqQual17  - SeqQual_sideband
   *
   *  for the called base, take the ratio of the value of that
   *  trace 1/2 of the way between this base and the next base
   *  over the value of the trace at its peak...
   *  compare that to the ratio of the value of that trace 1/2 or
   *  the way between this base and the previous base
   *  over the value of the trace at its peak...take the worst ratio
   *  and average that over the 8 bases before and after this base
   *
   */

#ifdef ANSI_C
static void SeqQual_sideband( PhredData *phredData )
#else
static void SeqQual_sideband( phredData )
PhredData *phredData;
#endif
{
  int i;
  int pos;
  int one_half_for;
  int one_half_back;
  int half_window_size = 8; /* 8 */
  int start_sum, end_sum;
  int j;
  FLOAT forward, backward;

  for( i = 0; i < phredData->numBase[FIL]; i++ )
  {
    phredData->qualIndex[i] = 0.0;
  }
  phredData->qualType = 17;


  for( j = 0; j < phredData->numBase[FIL]; j++ )
  {

    end_sum = j + half_window_size;
    start_sum = j - half_window_size;


    if( end_sum >= phredData->numBase[FIL] )
    {
      end_sum = phredData->numBase[FIL] - 2;
    }
    if( start_sum < 0 )
    {
      start_sum = 1;
    }


   for( i = start_sum; i <= end_sum; i++ )
   {

     one_half_for = one_half_forwards( phredData, i );
     one_half_back = one_half_backwards( phredData, i );

    pos = phredData->baseLoc[FIL][i];

    switch( phredData->base[FIL][i] )
    {
      case 'A':
      case 'a':
        if( phredData->trace[0][pos] > 0.0 )
        {
          forward  = phredData->trace[0][one_half_for] / phredData->trace[0][pos];
          backward = phredData->trace[0][one_half_back] / phredData->trace[0][pos];
          if( forward > backward )
          {
            phredData->qualIndex[j] += forward;
          }
          else
          {
            phredData->qualIndex[j] += backward;
          }
        }
        break;
      case 'C':
      case 'c':
        if( phredData->trace[1][pos] > 0.0 )
        {
          forward  = phredData->trace[1][one_half_for] / phredData->trace[1][pos];
          backward = phredData->trace[1][one_half_back] / phredData->trace[1][pos];
          if( forward > backward )
          {
            phredData->qualIndex[j] += forward;
          }
          else
          {
            phredData->qualIndex[j] += backward;
          }
        }
        break;
      case 'G':
      case 'g':
        if( phredData->trace[2][pos] > 0.0 )
        {
          forward  = phredData->trace[2][one_half_for] / phredData->trace[2][pos];
          backward = phredData->trace[2][one_half_back] / phredData->trace[2][pos];
          if( forward > backward )
          {
            phredData->qualIndex[j] += forward;
          }
          else
          {
            phredData->qualIndex[j] += backward;
          }
        }
        break;
      case 'T':
      case 't':
        if( phredData->trace[3][pos] > 0.0 )
        {
          forward  = phredData->trace[3][one_half_for] / phredData->trace[3][pos];
          backward = phredData->trace[3][one_half_back] / phredData->trace[3][pos];
          if( forward > backward )
          {
            phredData->qualIndex[j] += forward;
          }
          else
          {
            phredData->qualIndex[j] += backward;
          }
        }
        break;
      default:
        phredData->qualIndex[j] += 1.0;
        break;
    }
  }
      phredData->qualIndex[j] = phredData->qualIndex[j] / ( end_sum - start_sum + 1 );
  }

   /*
   ** identify quality index as type 17
   */
   phredData->qualType = 17;
}



#ifdef ANSI_C
static void findLeftQualCutoff( PhredData *phredData, int *start_point )
#else
static void findLeftQualCutoff( phredData, start_point )
PhredData *phredData;
int *start_point;
#endif
{
  int j;

  /*
  ** go calculate the quality measure ... this will take care
  ** of calculating it both for the the left and the right
  */

  /*
  ** had already called this for the overall trace quality
  */
  SeqQual_nonCalledOverCalled( phredData );

  /*
  ** start looking at start_point+STEP_SIZE because you really do not
  ** want to look at what comes before the cutoff.  You only care
  ** about what is after the cutoff...so in essence you are making
  ** your window centered on the left foot of the rectangle window rather
  ** than the center
  */

  for( j = *start_point + STEP_SIZE; j < phredData->numBase[FIL]; j += STEP_SIZE )
  {
/*
printf("phredData->qualIndex[%d] is %4.3f\n", j, phredData->qualIndex[j] );
*/
    if( phredData->qualIndex[j] < NONCALLED_OVER_CALLED_CUTOFF )
    {
      *start_point = j - STEP_SIZE;
      return;
    }
  }

  *start_point = j;
  return;

}



#ifdef ANSI_C
static int findRightQualCutoff( PhredData *phredData, int num_bases )
#else
static int findRightQualCutoff( phredData, num_bases )
PhredData *phredData;
int num_bases;
#endif
{
  int i, j;
  int num_problems; /* count of number of bases in this run having
                       a quality index above the cutoff */
  int num_allowed_problems = 4;  /* number of values above the
                                  cutoff allowed before setting
                                  the right cutoff */
  FLOAT cutoff = SIDEBAND_CUTOFF;
                 /* cutoff for the value of this side-band-ratio */
  int rightCutoff;
  int first_problem_base;


  /*
  ** go calculate the quality measure
  */
  SeqQual_sideband( phredData );

 
  /*
  ** step through all of the bases in STEP_SIZE increments
  */
  num_problems = 0;
  for( j = phredData->leftTrimPoint; j < phredData->numBase[FIL]; j += STEP_SIZE )
  {
    /*
    ** if the quality index exceeds the cutoff ... count it as a problem area
    */
    if( phredData->qualIndex[j] > cutoff )
    {
      num_problems++;
      if( num_problems == 1 )
      {
        first_problem_base = j;
      }
    }
    else
    {
      num_problems = 0;
    }

    /*
    ** if we have reached the num_allowed_problems over
    ** consecutive bases...then go ahead and assign
    ** the right cutoff to the point where the start
    ** of the problem bases was found
    */
    if( num_problems == num_allowed_problems )
    {
      rightCutoff = j - num_problems * STEP_SIZE;
    }
  }

  if( num_problems < num_allowed_problems )
  {
    rightCutoff = LAST_ALLOWED_BASE;
    if( rightCutoff >= num_bases )
    {
      rightCutoff = num_bases - 1;
    }
  }

  /*
  ** now go check the other quality measure from this cutoff
  ** base backwards....checking that we are not exceeding
  ** that rule...or could just take the more conservative
  ** of the two estimates....except then we're sometimes
  ** actually finding the left cutoff
  */
  SeqQual_nonCalledOverCalled( phredData );

  for( i = rightCutoff; i > phredData->leftTrimPoint; i -= 8 )
  {
    /*
    ** If two consecutive regions are good using the noncalled over called
    ** cutoff, then go ahead and set the right cutoff there
    */
    if( phredData->qualIndex[i] < NONCALLED_OVER_CALLED_CUTOFF && phredData->qualIndex[i-STEP_SIZE] < NONCALLED_OVER_CALLED_CUTOFF )
    {
      break;
    }
  }
  rightCutoff = i;
       


  /*
  ** ABSOLUTE CUTOFF IS LAST_ALLOWED_BASE
  */
  if( rightCutoff > LAST_ALLOWED_BASE )
  {
    rightCutoff = LAST_ALLOWED_BASE;
  }

  /* remember that the right cutoff in the trace structure is not
  ** the base position of the right cutoff...rather it is num_bases
  ** minus that position
  */
  rightCutoff = num_bases - rightCutoff;

  return( rightCutoff );
}



/*
** returns a one if the overall trace quality was good enough to
** warrant keeping the trace....0 if the trace should be thrown away

** Seq is the trace structure
** for step size...you are looking at values for the quality
** measure each STEP_SIZE-th base
*/
#ifdef ANSI_C
static int overallTraceQual( PhredData *phredData )
#else
static int overallTraceQual( phredData )
PhredData *phredData;
#endif

{
  int j;
  int num_good = 0; /* count of number of consecutive bases having
                     a quality index value better than the cutoff */
  int num_consecutive_good = 10;  /* 24 being about 200 bases of
                                   good trace since values are
                                   read every 8 bases */
  int num_problems = 0; /* count of number of bases in this run having
                       a quality index above the cutoff */
  int num_problems_allowed = 4;  /* number of values above the
                                  cutoff allowed in the span
                                  of num_consecutive_good */
  FLOAT cutoff = OVERALL_TRACE_QUAL_CUTOFF;  /* cutoff for the value of max_non_ca
lled over called */
  int last_problem;

  SeqQual_nonCalledOverCalled( phredData );

  for( j = 0; j < phredData->numBase[FIL]; j += STEP_SIZE )
  {
    if( phredData->qualIndex[j] < cutoff )
    {
      num_good++;
    }
    else
    {
      num_problems++;
      if( num_problems == 1 )
      {
        last_problem = j;
      }
      if( num_problems > num_problems_allowed )
      {
        num_good = 0;
        num_problems = 0;
        j = last_problem + 1;
      }
    }

    /*
    ** make sure you  hit num_consecutive_good in a row and that you are not
    ** out past LAST_ALLOWED_BASE when you do hit it
    */
    if( num_good == num_consecutive_good )
    {
      return( 1 );
    }
  }

  return( 0 );
}



/*
** looks for left cutoff, if it doesn't find a "enzInString", then
** it looks from enzInString less it's last character, etc 
*/
#ifdef ANSI_C
static int findLeftCutoff( PhredData *phredData, char *enzInString )
#else
static int findLeftCutoff( phredData, enzInString )
PhredData *phredData;
char *enzInString;
#endif
{
  int maxStartPos = 100; /* if the enzyme site wasn't found before this
                          baseNum, then that's probably not the cloning
                          site */
  int i, j, found;
  int indices[100];
  int num_matches;
  char *theSeq;
  int num_bases;
  char enzString[PHRED_MAX_STRING_LEN];
  int cut_point;

  cut_point = 0;
  found = 0;


  /*
  ** mdh - init array to known state (watch out for hardcoded 100).
  */
  memset( (void *)indices, 0, (size_t) ( sizeof(int)*100 ) );

  pstrcpy( enzString, enzInString, PHRED_MAX_STRING_LEN );

  num_bases = phredData->numBase[FIL];
  theSeq = (char *)calloc( num_bases, sizeof( char ) );

  j = 0;
  for( i = 0; i < num_bases; i++ )
  {
    theSeq[i] = phredData->base[FIL][i];
  }

  /*
  ** look for first occurrence of enzString;
  ** just look a match with at most i mismatches, starting
  ** with 0 mismatches down to two
  */
  for( i = 0; i < 3; i++ )
  {
    num_matches = stringMatch( enzString, strlen( enzString ), theSeq, num_bases, i, indices );
    if( num_matches > 0 )
    {
      if( indices[0] < maxStartPos )
      {
        found = 1;
        break;
      }
    }
  }

  free( theSeq );

  if( found && indices[0] )
  {
    cut_point = indices[0] + strlen( enzString );
  }
  else
  {
    cut_point = 0;

    /*
    ** make sure there are not a bunch of Ns from the ABI
    ** primer problem at the start of this sequence...move the
    ** left cutoff past all of those Ns
    */
    for( i = 0; i < phredData->numBase[FIL]; i++ )
    {
      if( phredData->base[FIL][i] != 'N' && phredData->base[FIL][i] != '-' )
      {
        cut_point = i;
        break;
      }
    }
  }

  findLeftQualCutoff( phredData, &cut_point );

  return( cut_point );
}



#ifdef ANSI_C
static int findRightCutoff( PhredData *phredData )
#else
static int findRightCutoff( phredData )
PhredData *phredData;
#endif
{
  int num_bases;
  char *theSeq;
  int i, j;
  int rightCutoff;

  num_bases = phredData->numBase[FIL];
  theSeq = (char *)calloc( num_bases, sizeof( char ) );

  j = 0;
  for( i = 0; i < num_bases; i++ )
  {
    theSeq[i] = phredData->base[FIL][i];
  }

  rightCutoff = findRightQualCutoff( phredData, num_bases );

  /*
  ** added so that the left and right cutoffs do not overlap
  */
  if( rightCutoff > num_bases - phredData->leftTrimPoint )
  {
     rightCutoff = num_bases - phredData->leftTrimPoint;
  }

  free( theSeq );

  return( rightCutoff );
}

/*
** trim sequence
*/
#ifdef ANSI_C
int trimSeq( PhredData *phredData, Option *option )
#else
int trimSeq( phredData, option )
PhredData *phredData;
Option *option;
#endif
{
  char *enzString;

  enzString = option->enzName;

  /*
  ** check the overall trace quality...it will
  ** return a zero if the trace should be thrown away
  */
  if( overallTraceQual( phredData ) )
  {
    phredData->leftTrimPoint = findLeftCutoff( phredData, enzString );
    phredData->rghtTrimPoint = findRightCutoff( phredData );
  }
  else
  {
    phredData->rghtTrimPoint = phredData->numBase[FIL];
    phredData->leftTrimPoint = 0;
  }

  return( OK );
}



