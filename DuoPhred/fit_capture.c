/** fit_capture.c **/

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
**    * fit_capture functions are based on ideas of Mike Wendl at GSC.  *    **
**                                                                           **
*******************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "phred.h"


/******************************************************************************/
/**                         CAPTURE_MISSED_PEAKS                             **/
/******************************************************************************/
#ifdef ANSI_C
int capture_missed_peaks1( int tr_length, Peak *first_peak,
                          Observed_peak *first_obs_peak, FLOAT **sig_func,
                          int begGoodTrace, int endGoodTrace )
#else
int capture_missed_peaks1( tr_length, first_peak, first_obs_peak, sig_func,
                           begGoodTrace, endGoodTrace )
int tr_length;
Peak *first_peak;
Observed_peak *first_obs_peak;
FLOAT **sig_func;
int begGoodTrace;
int endGoodTrace;
#endif
{
  int i, j;
  int i_called;
  int j_called, j_uncalled, j_most;
  FLOAT pred_period;
  FLOAT location;
  FLOAT frac_disp;
  FLOAT prev_period;
  FLOAT this_period;
  FLOAT next_period;
  Peak *peak, *captured_peak, *current_peak;

  /*
  ** SIFT THROUGH LIST OF BEST_UNCALLED_PEAKS TO SEE WHICH PEAKS WERE MISSED
  */
  for( peak = first_peak->next; peak; peak = peak->next )
  {
    if( peak->best_uncalled_peak && peak->obs_peak )
    {
      i = peak->best_uncalled_peak->i_maxx;
      j_uncalled = peak->best_uncalled_peak->nuc;
      i_called = peak->obs_peak->i_maxx;
      j_called = peak->obs_peak->nuc;
      
      if( !( sig_func[j_uncalled][i] >= 2.0 * sig_func[j_called][i_called] &&
          sig_func[j_called][i_called] < sig_func[j_uncalled][i_called] ) )
      {
        /*
        ** Crest of uncalled peak must be above other three traces.
        */
        j_most = 0;
        for( j = 1; j < 4; j++ )
        {
          if( sig_func[j][i] >= sig_func[j_most][i] )
          {
            j_most = j;
          }
        }
        if( j_uncalled != j_most )
        {
           continue;
        }
        else
        {
          j = j_uncalled;
        }				

        /*
        ** Must be a large peak.
        */
        if( sig_func[j_uncalled][i] < 0.33 * MAX_TRACE_VAL )
        {
          continue;
        }

        /*
        ** DO NOT BREAK UP SPLIT PEAKS
        */
        if( peak->obs_peak->type != 1 )
        {
          continue;
        }

        /*
        ** CHECK SPACING BETWEEN ADJACENT PEAKS
        */
        if( peak->best_uncalled_peak->i_maxx > peak->obs_peak->i_maxx )
        {
          if( peak->next && peak->next->obs_peak && peak->next->next && peak->next->next->obs_peak )
          {
            next_period = peak->next->next->obs_peak->split[1][1] -
                          peak->next->obs_peak->split[1][1];
          }
          else
          {
            next_period = 0.0;
          }
          if( peak->prev && peak->prev->obs_peak )
          {
            prev_period = peak->obs_peak->split[1][1] -
                            peak->prev->obs_peak->split[1][1];
          }
          else
          {
            prev_period = 0.0;
          }
          if( !peak->obs_peak || !peak->next || !peak->next->obs_peak )
          {
            continue;
          }
          if( peak->obs_peak->split[1][1] == peak->next->obs_peak->split[1][1] )
          {
            continue;
          }
          this_period = peak->next->obs_peak->split[1][1] -
                          peak->obs_peak->split[1][1];
          frac_disp = fabs( ( ( peak->obs_peak->split[1][1] +
                                 peak->next->obs_peak->split[1][1] ) / 2.0 -
                              peak->best_uncalled_peak->split[1][1] ) /
                            ( ( peak->obs_peak->split[1][1] -
                                peak->next->obs_peak->split[1][1] ) / 2.0 ) );

          if( frac_disp > 0.25 )
          {
            continue;
          }

          if( ( next_period == 0 && prev_period == 0 ) ||
              ( next_period > 0.0 &&
                ( next_period  - this_period * 0.5 ) / next_period > 0.40 ) &&
              ( prev_period > 0.0 &&
                ( prev_period  - this_period * 0.5 ) / prev_period > 0.40 ) )
          {
            continue;
          }

          if( peak->best_uncalled_peak->relative_area <
              peak->obs_peak->relative_area * 0.50 &&
              peak->best_uncalled_peak->relative_area <
              peak->next->obs_peak->relative_area * 0.50 )
          {
            continue;
          }
        }
        else
        {
          if( peak->next && peak->next->obs_peak )
          {
            next_period = peak->next->obs_peak->split[1][1] -
                          peak->obs_peak->split[1][1];
          }
          else
          {
            next_period = 0.0;
          }
          if( peak->prev && peak->prev->obs_peak && peak->prev->prev && peak->prev->prev->obs_peak )
          {
            prev_period = peak->prev->obs_peak->split[1][1] -
                          peak->prev->prev->obs_peak->split[1][1];
          }
          else
          {
            prev_period = 0.0;
          }
          if( !peak->obs_peak || !peak->prev || !peak->prev->obs_peak )
          {
            continue;
          }

          if( peak->obs_peak->split[1][1] == peak->prev->obs_peak->split[1][1] )
          {
            continue;
          }
          this_period = peak->obs_peak->split[1][1] -
                        peak->prev->obs_peak->split[1][1];
          frac_disp = fabs( ( ( peak->obs_peak->split[1][1] +
                                peak->prev->obs_peak->split[1][1] ) / 2.0 -
                              peak->best_uncalled_peak->split[1][1] ) /
                            ( ( peak->obs_peak->split[1][1] -
                                peak->prev->obs_peak->split[1][1] ) / 2.0 ) );
          if( frac_disp > 0.25 )
          {
            continue;
          }

          if( ( next_period == 0 && prev_period == 0 ) ||
              ( next_period > 0.0 &&
                ( next_period  - this_period * 0.5 ) / next_period > 0.40 ) &&
              ( prev_period > 0.0 &&
                ( prev_period  - this_period * 0.5 ) / prev_period > 0.40 ) )

          {
            continue;
          }

          if( peak->best_uncalled_peak->relative_area <
              peak->obs_peak->relative_area * 0.50 &&
              peak->best_uncalled_peak->relative_area <
              peak->prev->obs_peak->relative_area * 0.50 )
          {
            continue;
          }
        }

        if( i < begGoodTrace || i > 0.75 * endGoodTrace )
        {
          continue;
        }

        /*
        ** INSERT THIS PEAK INTO LIST TO BE CALLED SINCE IT PASSED FILTERS
        */
        captured_peak = alloc_predicted_peak();
        location = (FLOAT)peak->best_uncalled_peak->i_maxx;
        pred_period = peak->pred_period;
        init_peak(captured_peak, location, pred_period, 0.0, 0.0 );
        captured_peak->obs_peak = peak->best_uncalled_peak;
        captured_peak->obs_peak->peak[1] = captured_peak;
        captured_peak->obs_peak->type    = 1;
        peak->best_uncalled_peak = 0;
        if( captured_peak->obs_peak->i_maxx > peak->obs_peak->i_maxx )
        {
          current_peak = peak;
        }
        else
        {
          current_peak = peak->prev;
        }
        captured_peak->prev = current_peak;
        captured_peak->next = current_peak->next;
        if (current_peak->next) current_peak->next->prev = captured_peak;
        current_peak->next = captured_peak;
        peak = captured_peak;
      }
    }
  }

  return( OK );
}

typedef struct
{
  int   pos;
  int   type;
  int   base;
  int   capture;
  int   basMaxSpc;
  int   scrn;
  FLOAT relArea;
  FLOAT space;
  FLOAT spcRatio;
  FLOAT maxDownF1;
  FLOAT maxDown;
  Peak  *peak;
} EvalTrace;

/*
** Calculate space ratio in specified interval.
*/
#ifdef ANSI_C
static int calcSpcRatio( int numBase, EvalTrace *evalTrace, int beg, int end )
#else
static int calcSpcRatio( numBase, evalTrace, beg, end )
int numBase;
EvalTrace *evalTrace;
int beg;
int end;
#endif
{
  int i, j;
  int ibeg, iend;
  FLOAT min, max;

  /*
  ** Error on impossible request.
  */
  if( end < 3 || beg > numBase - 3 )
  {
    fprintf( stderr, "calcSpcRatio: error: invalid values: numBase: %d  begin: %d  end: %d\n",
                     numBase, beg, end );
    return( ERROR );
  }

  /*
  ** Adjust start and finish as required.
  */
  ibeg = ( beg >= 3 ) ? beg : 3;
  iend = ( end <= numBase - 3 ) ? end : numBase - 3;

  /*
  ** Calculate space ratio.
  */
  for( i = ibeg; i <= iend; ++i )
  {
    min = evalTrace[i-3].space;
    max = evalTrace[i-3].space;
    evalTrace[i].basMaxSpc = i - 3;
    for( j = -2; j <= 2; ++j )
    {
      if( evalTrace[i+j].space < min )
      {
        min = evalTrace[i+j].space;
      }
      if( evalTrace[i+j].space > max )
      {
        max = evalTrace[i+j].space;
        evalTrace[i].basMaxSpc = i + j;
      }
    }
    evalTrace[i].spcRatio = min > 0.0 ? max / min : 100.0;
  }

  if( beg < 3 )
  {
    for( i = beg; i < 3; ++i )
    {
      evalTrace[i].spcRatio  = evalTrace[3].spcRatio;
      evalTrace[i].basMaxSpc = evalTrace[3].basMaxSpc;
    }
  }
      
  if( end > numBase - 3 )
  {
    for( i = numBase - 2; i <= end; ++i )
    {
      evalTrace[i].spcRatio  = evalTrace[numBase-3].spcRatio;
      evalTrace[i].basMaxSpc = evalTrace[numBase-3].basMaxSpc;
    }
  }

  return( OK );
}

#ifdef ANSI_C
static FLOAT recalcSpcRatio( int numBase, int i, int k,
                      EvalTrace *evalTrace, int pos )
#else
static FLOAT recalcSpcRatio( numBase, i, k, evalTrace, pos )
int numBase;
int i;
int k;
EvalTrace *evalTrace;
int pos;
#endif
{
  int j;
  int beg, end;
  FLOAT spc1, spc2;
  FLOAT min, max;
  FLOAT spcRatio;

  beg = ( i - 3 >= 0 ) ? i - 3 : 0;
  end = ( i + 2 <= numBase - 1 ) ? i + 2 : numBase - 1;

  min = ( beg == k ) ? evalTrace[beg+1].space : evalTrace[beg].space;
  max = min;
  for( j = beg; j <= end; ++j )
  {
    if( j == k )
    {
      spc1 = (FLOAT)( pos - evalTrace[j].pos );
      spc2 = (FLOAT)( evalTrace[j+1].pos - pos );
      if( spc1 < min )
      {
        min = spc1;
      }
      if( spc1 > max )
      {
        max = spc1;
      }
      if( spc2 < min )
      {
        min = spc2;
      }
      if( spc2 > max )
      {
        max = spc2;
      }
    }
    else
    {
      if( evalTrace[j].space < min )
      {
        min = evalTrace[j].space;
      }
      if( evalTrace[j].space > max )
      {
        max = evalTrace[j].space;
      }
    }
  }

  spcRatio = ( min != 0.0 ) ? max / min : 100.0;

  return( spcRatio );
}


#ifdef ANSI_C
static Observed_peak *find_best_obs_peak( int numBase, EvalTrace *evalTrace,
                                   int i, int k, Observed_peak *first_obs_peak,
                                   int lo, int hi, FLOAT *spcRatio )
#else
static Observed_peak *find_best_obs_peak( numBase, evalTrace, i, k, first_obs_peak,
                                   lo, hi, spcRatio )
int numBase;
EvalTrace *evalTrace;
int i;
int k;
Observed_peak *first_obs_peak;
int lo;
int hi;
FLOAT *spcRatio;
#endif
{
  FLOAT rat;
  FLOAT minRat;
  Observed_peak *obs_peak;
  Observed_peak *best_obs_peak;

  minRat = evalTrace[i].spcRatio;
  best_obs_peak = NULL;

  /*
  ** Loop through observed peaks.
  */
  for( obs_peak = first_obs_peak; obs_peak && obs_peak->split[1][1] < hi; obs_peak = obs_peak->next )
  {
    /*
    ** Does peak meet candidate criteria?
    */
    if( obs_peak->split[1][1] > lo &&
        obs_peak->relative_area > 0.5 &&
        obs_peak->type == 0 )
    {
      rat = recalcSpcRatio( numBase, i, k, evalTrace, obs_peak->split[1][1] );
      if( rat < minRat )
      {
        best_obs_peak = obs_peak;
        minRat = rat;
      }
    }
  }
  *spcRatio = minRat;
  return( best_obs_peak );
}




/******************************************************************************/
/**                         CAPTURE_MISSED_PEAKS V2                          **/
/******************************************************************************/
#ifdef ANSI_C
int capture_missed_peaks2( int tr_length, FLOAT **tr_vals, FLOAT *tot_vals,
                          Peak *first_peak, Observed_peak *first_obs_peak )
#else
int capture_missed_peaks2( tr_length, tr_vals, tot_vals, first_peak, first_obs_peak )
int tr_length;
FLOAT **tr_vals;
FLOAT *tot_vals;
Peak *first_peak;
Observed_peak *first_obs_peak;
#endif
{
  int i, j, k;
  int itmp;
  int type;
  int nbyte;
  int shift;
  int numBase;
  int pos1, bas1;
  int pos2, bas2;
  FLOAT minPair;
  FLOAT location;
  FLOAT pred_period;
  FLOAT rat1, rat2;
  Peak *peak;
  Peak *current_peak;
  Peak *captured_peak;
  Observed_peak *obs_peak;
  EvalTrace *evalTrace;

  /*
  ** Count bases.
  */
  numBase = 0;
  for( peak = first_peak; peak; peak = peak->next )
  {
    ++numBase;
  }

  if( numBase < 10 )
  {
    return( OK );
  }

  /*
  ** Allocate memory.
  */
  nbyte = numBase * sizeof( EvalTrace );
  if( ( evalTrace = (EvalTrace *)malloc( nbyte ) ) == NULL )
  {
    fprintf( stderr, "unable to allocate memory\n" );
    return( ERROR );
  }

  /*
  ** Collect some information for trace evaluation.
  */
  i = 0;
  for( peak = first_peak; peak; peak = peak->next )
  {
    obs_peak = peak->obs_peak;
    if( obs_peak )
    {
      type = obs_peak->type;
      for( shift = 1; shift <= type; ++shift )
      {
        evalTrace[i].base    = obs_peak->nuc;
        evalTrace[i].pos     = obs_peak->split[type][shift];
        evalTrace[i].type    = type;
        evalTrace[i].relArea = obs_peak->relative_area;
        evalTrace[i].peak    = peak;
        evalTrace[i].capture = 0;
        ++i;
        if( shift != type )
        {
          peak = peak->next;
        }
      }
    }
    else
    {
      evalTrace[i].base = 4;
      evalTrace[i].pos  = peak->pred_location;
      evalTrace[i].type = 0;
      evalTrace[i].relArea = 0.0;
      evalTrace[i].peak = peak;
      evalTrace[i].capture = 0;
      ++i;
    }
  }

  numBase = i;

  /*
  ** Calculate space between pairs of bases.
  */
  for( i = 0; i < ( numBase - 1 ); ++i )
  {
    evalTrace[i].space = evalTrace[i+1].pos - evalTrace[i].pos;
  }
  evalTrace[numBase-1].space = evalTrace[numBase-2].space;

  /*
  ** Calculate the low points in the intervals between peaks.
  */
  for( i = 0; i < ( numBase - 1 ); ++i )
  {
    pos1 = evalTrace[i].pos < tr_length ? evalTrace[i].pos : tr_length - 1;
    bas1 = evalTrace[i].base;
    pos2 = evalTrace[i+1].pos < tr_length ? evalTrace[i+1].pos : tr_length - 1;
    bas2 = evalTrace[i+1].base;
    if( pos1 > pos2 )
    {
      itmp = pos1;
      pos1 = pos2;
      pos2 = itmp;

      itmp = bas1;
      bas1 = bas2;
      bas2 = itmp;
    }
    if( bas1 < 4 && bas2 < 4 )
    {
      minPair = ( tr_vals[bas1][pos1] < tr_vals[bas2][pos2] ) ?
                tr_vals[bas1][pos1] : tr_vals[bas2][pos2];
    }
    else
    if( bas1 == 4 && bas2 == 4 )
    {
      minPair = -1.0;
    }
    else
    if( bas1 == 4 )
    {
      minPair = tr_vals[bas2][pos2];
    }
    else
    if( bas2 == 4 )
    {
      minPair = tr_vals[bas1][pos1];
    }
    else
    {
      printf( "unrecognized base type\n" );
      minPair = 0.0;
    }

    k = pos1;
    for( j = pos1; j <= pos2; ++j )
    {
      if( tot_vals[j] < tot_vals[k] )
      {
        k = j;
      }
    }
    evalTrace[i].maxDownF1 = tot_vals[k];

    if( minPair > 0.0 )
    {
      evalTrace[i].maxDownF1 /= minPair;
    }
    else
    {
      evalTrace[i].maxDownF1 = -1.0;
    }
  }
  evalTrace[numBase-1].maxDownF1 = evalTrace[numBase-2].maxDownF1;

  /*
  ** Calculate the maxDown count.
  */
  for( i = 0; i < numBase; ++i )
  {
    j = 1;
    if( evalTrace[i].relArea > 0.0 &&
        evalTrace[i].type == 1 &&
        evalTrace[i].base != 'N' )
    {
      while( 1 )
      {
        if( i - j < 0 ||
            evalTrace[i-j].maxDownF1 > 0.99 ||
            evalTrace[i-j].relArea <= 0.0 ||
            evalTrace[i-j].type != 1 ||
            evalTrace[i-j].base == 'N' )
        {
          break;
        }
        if( i + j >= numBase ||
            evalTrace[i+j-1].maxDownF1 > 0.99 ||
            evalTrace[i+j].relArea <= 0.0 ||
            evalTrace[i+j].type != 1 ||
            evalTrace[i+j].base == 'N' )
        {
          break;
        }
        ++j;
      }
    }
    evalTrace[i].maxDown = (float)( -( j - 1 ) );
  }

  /*
  ** Calculate the space ratio.
  */
  if( calcSpcRatio( numBase, evalTrace, 0, numBase - 1 ) == ERROR )
  {
    fprintf( stderr, "capture_missed_peaks: bad status: calcSpcRatio\n" );
    free( evalTrace );
    return( ERROR );
  }

  /*
  ** Screen out split peaks in the window.
  */
  for( i = 3; i < numBase - 3; ++i )
  {
    evalTrace[i].scrn = 0;
    for( j = -3; j <= 3; ++j )
    {
      if( evalTrace[i+j].type != 1 )
      {
        evalTrace[i].scrn = 1;
      }
    }
  }
  evalTrace[2].scrn = evalTrace[3].scrn;
  evalTrace[1].scrn = evalTrace[3].scrn;
  evalTrace[0].scrn = evalTrace[3].scrn;
  evalTrace[numBase-3].scrn = evalTrace[numBase-4].scrn;
  evalTrace[numBase-2].scrn = evalTrace[numBase-4].scrn;
  evalTrace[numBase-1].scrn = evalTrace[numBase-4].scrn;

  /*
  ** Try to capture missed peaks.
  */
  for( i = 0; i < numBase - 1; ++i )
  {
    /*
    ** Is this group of peaks a candidate for peak capturing?
    */
    k = evalTrace[i].basMaxSpc;
    if( k == i &&
        evalTrace[i].maxDown <= -4.0 &&
        evalTrace[i].spcRatio >= 1.5 &&
        evalTrace[i].spcRatio <  3.0 &&
        evalTrace[i].scrn == 0 )
    {
      /*
      ** Examine the big space in the group of peaks.
      ** Is there a suitable uncalled peak?
      */
      peak = 0;
      current_peak = evalTrace[k].peak;

      if( evalTrace[k].capture == 0 &&
          evalTrace[k].peak->best_uncalled_peak &&
          evalTrace[k].peak->best_uncalled_peak->relative_area >= 0.5 &&
          evalTrace[k].peak->best_uncalled_peak->split[1][1] >
          evalTrace[k].peak->obs_peak->split[1][1] &&
          evalTrace[k].peak->best_uncalled_peak->split[1][1] <
          evalTrace[k+1].peak->obs_peak->split[1][1] &&

          evalTrace[k+1].peak->best_uncalled_peak &&
          evalTrace[k+1].peak->best_uncalled_peak->relative_area >= 0.5 &&
          evalTrace[k+1].peak->best_uncalled_peak->split[1][1] <
          evalTrace[k+1].peak->obs_peak->split[1][1] &&
          evalTrace[k+1].peak->best_uncalled_peak->split[1][1] >
          evalTrace[k].peak->obs_peak->split[1][1] )
      {
        /*
        ** Two suitable uncalled peaks.  Does peak capturing improve the
        ** space ratio?  Choose the peak with the best improvement.
        */
        rat1 = recalcSpcRatio( numBase, i, k, evalTrace,
                               evalTrace[k].peak->best_uncalled_peak->split[1][1] );
        rat2 = recalcSpcRatio( numBase, i, k, evalTrace,
                               evalTrace[k+1].peak->best_uncalled_peak->split[1][1] );
        if( rat1 < rat2 )
        {
          if( rat1 < 1.00 * evalTrace[i].spcRatio )
          {
            peak = evalTrace[k].peak;
          }
        }
        else
        {
          if( rat2 < 1.00 * evalTrace[i].spcRatio )
          {
            peak = evalTrace[k+1].peak;
          }
        }
      }
      else
      if( evalTrace[k].capture == 0 &&
          evalTrace[k].peak->best_uncalled_peak &&
          evalTrace[k].peak->best_uncalled_peak->relative_area >= 0.5 &&
          evalTrace[k].peak->best_uncalled_peak->split[1][1] >
          evalTrace[k].peak->obs_peak->split[1][1] &&
          evalTrace[k].peak->best_uncalled_peak->split[1][1] <
          evalTrace[k+1].peak->obs_peak->split[1][1] )
      {
        /*
        ** One suitable uncalled peak. Does peak capture improve the
        ** space ratio?
        */
        rat1 = recalcSpcRatio( numBase, i, k, evalTrace,
                               evalTrace[k].peak->best_uncalled_peak->split[1][1] );
        if( rat1 < 1.00 * evalTrace[i].spcRatio )
        {
          peak = evalTrace[k].peak;
        }
      }
      else
      if( evalTrace[k].capture == 0 &&
          evalTrace[k+1].peak->best_uncalled_peak &&
          evalTrace[k+1].peak->best_uncalled_peak->relative_area >= 0.5 &&
          evalTrace[k+1].peak->best_uncalled_peak->split[1][1] <
          evalTrace[k+1].peak->obs_peak->split[1][1] &&
          evalTrace[k+1].peak->best_uncalled_peak->split[1][1] >
          evalTrace[k].peak->obs_peak->split[1][1] )
      {
        /*
        ** One suitable uncalled peak. Does peak capture improve the
        ** space ratio?
        */
        rat2 = recalcSpcRatio( numBase, i, k, evalTrace,
                               evalTrace[k+1].peak->best_uncalled_peak->split[1][1] );
        if( rat2 < 1.00 * evalTrace[i].spcRatio )
        {
          peak = evalTrace[k+1].peak;
        }
      }

      if( peak == 0 )
      {
        continue;
      }

      /*
      ** Insert peak.
      */
      captured_peak = alloc_predicted_peak();
      location      = (FLOAT)peak->best_uncalled_peak->split[1][1];
      pred_period   = peak->pred_period;
      init_peak( captured_peak, location, pred_period, 0.0, 0.0 );

      captured_peak->obs_peak          = peak->best_uncalled_peak;
      captured_peak->obs_peak->peak[1] = captured_peak;
      captured_peak->obs_peak->type    = 1;

      captured_peak->next              = current_peak->next;
      current_peak->next               = captured_peak;
      captured_peak->prev              = current_peak;
      if( captured_peak->next )
      {
        captured_peak->next->prev      = captured_peak;
      }

      peak->best_uncalled_peak         = 0;

      evalTrace[k].capture = 1;

    }
  }

  /*
  ** Free memory.
  */
  free( evalTrace );

  return( OK );
}

#ifdef ANSI_C
int capture_missed_peaks3( int tr_length, FLOAT **tr_vals, FLOAT *tot_vals,
                           Peak *first_peak, Observed_peak *first_obs_peak,
                           int npk, LocPeak *locPeak )
#else
int capture_missed_peaks3( tr_length, tr_vals, tot_vals, first_peak, first_obs_peak,
                           npk, locPeak )
int tr_length;
FLOAT **tr_vals;
FLOAT *tot_vals;
Peak *first_peak;
Observed_peak *first_obs_peak;
int npk;
LocPeak *locPeak;
#endif
{
  int i, j, k;
  int itmp;
  int pos1, bas1;
  int pos2, bas2;
  int type;
  int nbyte;
  int shift;
  int numBase;
  int loc, lo, hi;
  int *pkMatch;
  FLOAT minPair;
  FLOAT location;
  FLOAT pred_period;
  Peak *peak;
  Peak *current_peak;
  Peak *captured_peak;
  Observed_peak *obs_peak;
  EvalTrace *evalTrace;

  /*
  ** Count bases.
  */
  numBase = 0;
  for( peak = first_peak; peak; peak = peak->next )
  {
    ++numBase;
  }

  /*
  ** Allocate memory.
  */
  nbyte = numBase * sizeof( EvalTrace );
  if( ( evalTrace = (EvalTrace *)malloc( nbyte ) ) == NULL )
  {
    fprintf( stderr, "unable to allocate memory\n" );
    return( ERROR );
  }

  nbyte = npk * sizeof( int );
  if( ( pkMatch = (int *)malloc( nbyte ) ) == NULL )
  {
    fprintf( stderr, "unable to allocate memory\n" );
    return( ERROR );
  }

  /*
  ** Collect some information for trace evaluation.
  */
  i = 0;
  for( peak = first_peak; peak; peak = peak->next )
  {
    obs_peak = peak->obs_peak;
    if( obs_peak )
    {
      type = obs_peak->type;
      for( shift = 1; shift <= type; ++shift )
      {
        evalTrace[i].base    = obs_peak->nuc;
        evalTrace[i].pos     = obs_peak->split[type][shift];
        evalTrace[i].type    = type;
        evalTrace[i].relArea = obs_peak->relative_area;
        evalTrace[i].peak    = peak;
        evalTrace[i].capture = 0;
        ++i;
        if( shift != type )
        {
          peak = peak->next;
        }
      }
    }
    else
    {
      evalTrace[i].base = 4;
      evalTrace[i].pos  = peak->pred_location;
      evalTrace[i].type = 0;
      evalTrace[i].relArea = 0.0;
      evalTrace[i].peak = peak;
      evalTrace[i].capture = 0;
      ++i;
    }
  }

  numBase = i;

  /*
  ** Calculate the low points in the intervals between peaks.
  */
  for( i = 0; i < ( numBase - 1 ); ++i )
  {
    pos1 = evalTrace[i].pos < tr_length ? evalTrace[i].pos : tr_length - 1;
    bas1 = evalTrace[i].base;
    pos2 = evalTrace[i+1].pos < tr_length ? evalTrace[i+1].pos : tr_length - 1;
    bas2 = evalTrace[i+1].base;
    if( pos1 > pos2 )
    {
      itmp = pos1;
      pos1 = pos2;
      pos2 = itmp;

      itmp = bas1;
      bas1 = bas2;
      bas2 = itmp;
    }
    if( bas1 < 4 && bas2 < 4 )
    {
      minPair = ( tr_vals[bas1][pos1] < tr_vals[bas2][pos2] ) ?
                tr_vals[bas1][pos1] : tr_vals[bas2][pos2];
    }
    else
    if( bas1 == 4 && bas2 == 4 )
    {
      minPair = -1.0;
    }
    else
    if( bas1 == 4 )
    {
      minPair = tr_vals[bas2][pos2];
    }
    else
    if( bas2 == 4 )
    {
      minPair = tr_vals[bas1][pos1];
    }
    else
    {
      printf( "unrecognized base type\n" );
      minPair = 0.0;
    }

    k = pos1;
    for( j = pos1; j <= pos2; ++j )
    {
      if( tot_vals[j] < tot_vals[k] )
      {
        k = j;
      }
    }
    evalTrace[i].maxDownF1 = tot_vals[k];

    if( minPair > 0.0 )
    {
      evalTrace[i].maxDownF1 /= minPair;
    }
    else
    {
      evalTrace[i].maxDownF1 = -1.0;
    }
  }
  evalTrace[numBase-1].maxDownF1 = evalTrace[numBase-2].maxDownF1;

  /*
  ** Calculate the maxDown count.
  */
  for( i = 0; i < numBase; ++i )
  {
    j = 1;
    if( evalTrace[i].relArea > 0.0 &&
        evalTrace[i].type == 1 &&
        evalTrace[i].base != 'N' )
    {
      while( 1 )
      {
        if( i - j < 0 ||
            evalTrace[i-j].maxDownF1 > 0.99 ||
            evalTrace[i-j].relArea <= 0.0 ||
            evalTrace[i-j].type != 1 ||
            evalTrace[i-j].base == 'N' )
        {
          break;
        }
        if( i + j >= numBase ||
            evalTrace[i+j-1].maxDownF1 > 0.99 ||
            evalTrace[i+j].relArea <= 0.0 ||
            evalTrace[i+j].type != 1 ||
            evalTrace[i+j].base == 'N' )
        {
          break;
        }
        ++j;
      }
    }
    evalTrace[i].maxDown = (float)( -( j - 1 ) );
  }

  /*
  ** Screen out split peaks in the window.
  */
  for( i = 3; i < numBase - 3; ++i )
  {
    evalTrace[i].scrn = 0;
    for( j = -3; j <= 3; ++j )
    {
      if( evalTrace[i+j].type != 1 )
      {
        evalTrace[i].scrn = 1;
      }
    }
  }
  evalTrace[2].scrn = evalTrace[3].scrn;
  evalTrace[1].scrn = evalTrace[3].scrn;
  evalTrace[0].scrn = evalTrace[3].scrn;
  evalTrace[numBase-3].scrn = evalTrace[numBase-4].scrn;
  evalTrace[numBase-2].scrn = evalTrace[numBase-4].scrn;
  evalTrace[numBase-1].scrn = evalTrace[numBase-4].scrn;


  /*
  ** Match evalTrace peaks to locPeak peaks.
  */
  for( i = 0; i < npk; ++i )
  {
    pkMatch[i] = -1;
  }
  for( i = 0; i < npk; ++i )
  {
    loc = locPeak[i].point;
    lo = ( loc - 3 ) >= 0 ? loc - 3 : 0;
    hi = ( loc + 3 ) < tr_length ? loc + 3 : tr_length - 1;
    for( j = 0; j < numBase; ++j )
    {
      if( evalTrace[j].pos >= lo && evalTrace[j].pos <= hi )
      {
        pkMatch[i] = j;
        break;
      }
    }
  }

  /*
  ** Capture peaks.
  */
  for( i = 1; i < npk - 1; ++i )
  {
    j = pkMatch[i];
    if( j == -1 )
    {
      if( locPeak[i].spcRatio < 2.0 &&
          pkMatch[i-1] > 0 &&
          pkMatch[i+1] > 0 &&
          evalTrace[pkMatch[i-1]].maxDown <= -3.0 &&
          evalTrace[pkMatch[i-1]].scrn == 0 &&
          evalTrace[pkMatch[i+1]].maxDown <= -3.0 &&
          evalTrace[pkMatch[i+1]].scrn == 0 )
      {
        /*
        ** Missing high quality peak.  Find corresponding observed peak and
        ** insert it.
        */
        loc = locPeak[i].point;
        lo = ( loc - 3 ) >= 0 ? loc - 3 : 0;
        hi = ( loc + 3 ) < tr_length ? loc + 3 : tr_length - 1;
        for( obs_peak = first_obs_peak; obs_peak; obs_peak = obs_peak->next )
        {
          if( obs_peak->type == 0 &&
              obs_peak->split[1][1] >= lo &&
              obs_peak->split[1][1] <= hi &&
              obs_peak->nuc == locPeak[i].base &&
              obs_peak->relative_area > 0.5 )
          {
            k = pkMatch[i-1];
            current_peak = evalTrace[k].peak;
            captured_peak = alloc_predicted_peak();
            location      = (FLOAT)obs_peak->split[1][1];
            pred_period   = evalTrace[k].peak->pred_period;
            init_peak( captured_peak, location, pred_period, 0.0, 0.0 );
            captured_peak->obs_peak          = obs_peak;
            captured_peak->obs_peak->peak[1] = captured_peak;
            captured_peak->obs_peak->type    = 1;

            captured_peak->next              = current_peak->next;
            current_peak->next               = captured_peak;
            captured_peak->prev              = current_peak;
            if( captured_peak->next )
            {
              captured_peak->next->prev        = captured_peak;
            }

            /*
            ** If obs_peak was a best_uncalled_peak, remove it.
            */
            for( peak = first_peak; peak; peak = peak->next )
            {
              if( peak->best_uncalled_peak == obs_peak )
              {
                peak->best_uncalled_peak = 0;
                break;
              }
            }

            break;
          }
        }
      }
    }
  }

  /*
  ** Free memory.
  */
  free( pkMatch );
  free( evalTrace );

  return( OK );
}

