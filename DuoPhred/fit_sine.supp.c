/** fit_sine.supp.c **/

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


/*
** These variables allow us to free dynamically
** allocated memory by calling free_predicted_peaks()
** and free_cos_sin().
*/
static int num_predicted_peaks;
static Peak *alloc_predicted_peaks;
static Cos_sin_array *alloc_cos_sin;

/*
** normalize tr_vals, so that all nucleotides have same
** average peak values (in a good part of the trace).
** The algorithm is not very robust (if there are lots
** of spurious peaks) and in any case may be totally
** unnecessary.
** tot_vals, the summed signal, is returned
*/
#ifdef ANSI_C
int normalize_traces( FILE *fp, FLOAT **tr_vals, FLOAT *tot_vals,
                      int tr_length, FLOAT *scale_traces )
#else
int normalize_traces( fp, tr_vals, tot_vals, tr_length, scale_traces )
FILE *fp;		/* pointer to .log file */
FLOAT **tr_vals;
FLOAT *tot_vals;
int tr_length;
FLOAT *scale_traces;
#endif
{
  int i, j, k;
  int lo, hi;
  int nuc_tr_counts[4];
  FLOAT ratio;
  FLOAT nuc_tr_sums[4];

  for( i = 0; i < tr_length; ++i )
  {
    tot_vals[i] = 0.0;
    for( j = 0; j < 4; j++ )
    {
      tot_vals[i] += tr_vals[j][i];
    }
  }

  for( j = 0; j < 4; j++ )
  {
    nuc_tr_counts[j] = 0;
    nuc_tr_sums[j] = 0.0;
  }

  /*
  ** Find peaks in total signal
  ** in central region of trace, and
  ** increment sum for appropriate nucleotide
  ** at each peak.
  */
  lo = 1500;
  hi = 2500;
  for( i = lo; i < hi; ++i )
  {
    if( tot_vals[i] > tot_vals[i-1] && tot_vals[i] > tot_vals[i+1] )
    {
      for( j = 0, k = 1; k < 4; k++ ) 
      {
	if( tr_vals[k][i] > tr_vals[j][i] )
        {
          j = k;
        }
      }
      nuc_tr_counts[j] += 1;
      nuc_tr_sums[j] += tr_vals[j][i];
    }
  }

  /*
  ** find average peak signal for each nucleotide
  */
  for( j = 0; j < 4; j++ )
  {
    if( nuc_tr_counts[j] )
    {
      nuc_tr_sums[j] /= nuc_tr_counts[j];
    }
    else
    {
      /*
      ** The polyData modules need useful scaling values.
      */
      scale_traces[j] = 0.0;
    }
    scale_traces[j] = nuc_tr_sums[j];
  }

  /*
  ** skip nucleotides which did not occur
  ** (possibly due to very weak signal)
  */

  /*
  ** find nucleotide with smallest average peak signal
  */
  for (j = -1, k = 0; k < 4; k++)
  {
    if( nuc_tr_counts[k] && ( j < 0 || nuc_tr_sums[k] < nuc_tr_sums[j] ) )
    {
      j = k;
    }
  }

  /*
  ** normalize signal to that of nucleotide
  ** with smallest average peak signal
  */
  for( k = 0; k < 4; k++ )
  { 
    if( nuc_tr_counts[k] )
    {
      ratio = nuc_tr_sums[j] / nuc_tr_sums[k];
      for( i = 0; i < tr_length; i++ )
      {
        tr_vals[k][i] = tr_vals[k][i] * ratio;
      }
    }
  }

  return( 0 );
}



#ifdef ANSI_C
int makeTotVals( FLOAT **tr_vals, FLOAT *tot_vals, int tr_length )
#else
int makeTotVals( tr_vals, tot_vals, tr_length )
FLOAT **tr_vals;
FLOAT *tot_vals;
int tr_length;
#endif
{

  int i, j;

  /*
  ** Find max.
  */
  for( i = 0; i < tr_length; i++ )
  {  
    tot_vals[i] = tr_vals[0][i];
    for( j = 1; j < 4; j++ )
    {
      if( tot_vals[i] < tr_vals[j][i] )
      {
        tot_vals[i] = tr_vals[j][i];
      }
    }
  }

  return( 0 );

}


/*
** make_predicted_peaks:
**
** Purpose
** =======
**
** Make a list of regularly spaced predicted peaks
** by fitting a sine curve to the summed trace values.
**
** Arguments
** =========
**
**         fp        (input)  FILE (*)
**                   Pointer to open log file.
**
**         tr_vals   (input)  FLOAT array, dimensions [4][tr_length]
**                   Trace values for the four bases.
**
**         tr_length (input)  int
**                   Number of points in trace.
**
**         tot_vals  (input)  FLOAT array, dimension [tr_length]
**                   Each value of tot_vals is the sum of the
**                   four trace values at that point.
**
**         make_predicted_peaks (output) linked list of Peak, dimension[]
**                              First peak in linked list of predicted peaks.
**
**
*/
/*
** bge
** added tr_vals to argument list
*/
#ifdef ANSI_C
Peak *make_predicted_peaks( FILE *fp, FLOAT **tr_vals, int tr_length, FLOAT *tot_vals,
                            int numLocPeak, LocPeak *locPeak,
                            int beginPeakPredictionOption,
                            int beginPeakPredictionPoint )
#else
Peak *make_predicted_peaks( fp, tr_vals, tr_length, tot_vals, numLocPeak, locPeak,
                            beginPeakPredictionOption, beginPeakPredictionPoint )
FILE *fp;
FLOAT **tr_vals;
int tr_length;
FLOAT *tot_vals;
int numLocPeak;
LocPeak *locPeak;
int beginPeakPredictionOption;
int beginPeakPredictionPoint;
#endif
{
  int i, k;
  int block_size, i_peak, min_i;
  int seg_length, flag;
  int site;
  FLOAT prev_ratio, prop_fitted;
  FLOAT best_offset, total_signal;
  FLOAT lower_period, upper_period;
  FLOAT old_best_period;
  FLOAT best_value, best_period, max_value;
  FLOAT middle_peak, peak_loc;
  FLOAT *temp_vals;
  Peak *predicted_peaks, *peak;
  Peak *first_peak, *last_nuc;

  int *tr_qual;
  FLOAT *syn_peak;
  FLOAT *ptr_sav;
  Option *option;

  option = getOption();

  /*
  **
  ** Store arrays used to restrict period variation and
  ** evaluate period.
  */
  tr_qual = (int *)ourMalloc( tr_length * sizeof( int ) );
  syn_peak = (FLOAT *)ourMalloc( tr_length * sizeof( FLOAT ) );

  /*
  ** initialize and allocate memory
  */
  block_size = 200;

  temp_vals = (FLOAT *)ourMalloc( tr_length * sizeof( FLOAT ) );
  predicted_peaks = (Peak *)ourMalloc( tr_length * sizeof( Peak ) );
  alloc_predicted_peaks = predicted_peaks;
  prev_ratio = 0.0;

  /*
  ** Find a good place to begin predicting peaks and
  ** generate "synthetic peaks", a square wave function
  ** that closely corresponds to the peak locations -
  ** this removes the effect of varying peak amplitude
  ** and creates a more symmetric function for the
  ** Fourier analysis.  Also create an evaluation of
  ** the trace that's used to prevent period variation
  ** in regions of poor trace quality.
  */
  i = findStartPred( tr_length, numLocPeak, locPeak, block_size, tr_qual );


  if( beginPeakPredictionOption )
  {
    if( beginPeakPredictionPoint > ( block_size / 2 )  &&
        beginPeakPredictionPoint < ( tr_length - block_size / 2 ) )
    {
      i = beginPeakPredictionPoint;
    }
    else
    {
      fprintf( stderr,
               "  invalid peak prediction start point: %d\n",
               beginPeakPredictionPoint );
      fprintf( stderr,
               "  must be between %d and %d: using default\n",
               block_size / 2 + 1,
               tr_length - 1 - block_size / 2 );
    }
  }

  if( option->verboseOption == 1 &&
      option->verboseLevel >= 63 )
  {
    fprintf( stderr,
             "make_predicted_peaks: begin peak prediction at point %d\n",
             i );
  }

  makePeakSignal( tr_length, numLocPeak, locPeak, syn_peak );

  ptr_sav = tot_vals;
  tot_vals = syn_peak;

  /*
  ** Apply a digital filter to the block of values.
  */
  total_signal = weight_vec( temp_vals, tot_vals + i - block_size / 2,
                             block_size, block_size / 2 );

  /*
  ** Estimate peak period using Fourier transform.
  */
  max_from_fft( temp_vals, block_size,
                &best_period, &best_value, 5.0, 20.0 );

  if( option->verboseOption == 1 &&
      option->verboseLevel >= 63 )
  {
    fprintf( stderr,
             "make_predicted_peaks: initial peak spacing is %f\n",
             best_period );
  }

  /*
  ** Sometimes 'findStartPred' returns a bad location (or the trace is
  ** lousy), which causes max_from_fft to find a ridiculous best_period.
  ** As a quick fix, just set the best_period to a nominal value in this
  ** case.
  */
  if( best_period > 50.0 )
  {
    best_period = 12.0;
  }

  /*
  ** proportion of the filtered peaks fitted by
  ** sinusoid of best_period.
  */
  prop_fitted = best_value / total_signal;
  best_period = (int)( 10.0 * best_period + .5 ) / 10.0;

  /*
  ** vary the period slightly and look for better
  ** periodicity and calculate the phase offset.
  */
  find_best( temp_vals, block_size,
             best_period, best_period, best_period, 0.1,
             &max_value, &best_period, &best_offset );
  
  /*
  ** find peak nearest to point i
  */
  best_offset += i - block_size / 2;
  middle_peak = nearest_peak( (FLOAT)i, best_offset, best_period );
  
  /*
  ** store peak information in predicted_peaks
  */
  init_peak( predicted_peaks, middle_peak, best_period,
             total_signal, prop_fitted );

  i_peak = 1;

  last_nuc = 0;
  best_period = predicted_peaks->pred_period;
  old_best_period = best_period;
  peak = predicted_peaks;

  /*
  ** Predict peaks, starting in middle
  ** (where good-fitting sine curve was found)
  ** and working toward back.
  */
  for( k = middle_peak + best_period + .5;
       k < tr_length; 
       k = peak_loc + best_period + .5 )
  {
    seg_length = block_size;
    if( k + seg_length / 2 <= tr_length )
    {
      min_i = k - seg_length / 2;
      site = seg_length / 2;
    }
    else
    {
      min_i = tr_length - seg_length;
      site = ( k < tr_length - 50 ) ? k - min_i : seg_length - 50;
    }
    total_signal = weight_vec( temp_vals, tot_vals + min_i,
                               seg_length, site );

    if( tr_qual[k] == 0 )
    {
      old_best_period = best_period;
      if( last_nuc )
      {
        lower_period = upper_period = best_period;
      }
      else
      {
        /*
        ** 1.5
        */
        lower_period = best_period - .05;
        upper_period = best_period + .05;
      }
    }
    else
    {

      lower_period = upper_period = best_period = old_best_period;

    }

    /*
    ** look for best-fitting sine curve,
    ** considering only periods which don't
    ** differ too much from the period at
    ** the previous base
    */
    flag = find_best( temp_vals, seg_length,
                      best_period, lower_period, upper_period,
                      .03, &max_value, &best_period, &best_offset );

    /*
    ** locate predicted peak nearest to k
    */
    best_offset += min_i;
    peak_loc = nearest_peak( (FLOAT)k, best_offset, best_period );
    if( peak_loc >= tr_length )
    {
      break;
    }
    peak = predicted_peaks + i_peak;

    /*
    ** store peak information
    */
    init_peak( peak, peak_loc, best_period,
               total_signal, max_value / total_signal );

    /*
    ** link peak into list
    */
    (peak - 1)->next = peak;

    ++i_peak;
  }

  /*
  ** Now start at middle and work toward front of trace.
  */
  best_period = predicted_peaks->pred_period;
  old_best_period = best_period;
  peak->next = 0;
  for( k = middle_peak - best_period + .5;
       k > 0;
       k = peak_loc - best_period + .5 )
  {
    seg_length = block_size;
    if( k - seg_length / 2 >= 0 )
    {
      min_i = k - seg_length / 2;
      site = seg_length / 2;
    }
    else
    {
      min_i = 0;
      site = ( k > 50 ) ? k : 50;
    }
    total_signal = weight_vec( temp_vals, tot_vals + min_i,
                               seg_length, site );

    if( tr_qual[k] == 0 )
    {
      old_best_period = best_period;
    }
    else
    {
      best_period = old_best_period;
    }

    flag = find_best( temp_vals, seg_length,
                      best_period, best_period - 0.25, best_period + 0.25,
                      0.05, &max_value, &best_period, &best_offset );
/*
    flag = find_best( temp_vals, seg_length,
                      best_period, best_period - 0.05, best_period + 0.05,
                      0.01, &max_value, &best_period, &best_offset );
*/

    best_offset += min_i;
    peak_loc = nearest_peak( (FLOAT)k, best_offset, best_period );
    if( peak_loc < 0.0 )
    {
      break;
    }
    peak = predicted_peaks + i_peak;
    init_peak( peak,  peak_loc, best_period,
               total_signal, max_value / total_signal );
    peak->next = (peak - 1)->next ? peak - 1 : predicted_peaks;
    i_peak++;
  }

  /*
  ** linked list details
  */
  first_peak = peak;
  first_peak->prev = 0;
  for (peak = first_peak; peak; peak = peak->next)
  {
    /*
    ** complete double linking
    */
    if( peak->next )
    {
      peak->next->prev = peak;
    }
  }

  num_predicted_peaks = i_peak;

  tot_vals = ptr_sav;
  ourFree( (char *)tr_qual );
  ourFree( (char *)syn_peak );

  ourFree( (char *)temp_vals );
  free_cos_sin();

  return( first_peak );
}  



#ifdef ANSI_C
Peak *alloc_predicted_peak( void )
#else
Peak *alloc_predicted_peak()
#endif
{
  Peak *next_predicted_peak;

  next_predicted_peak = alloc_predicted_peaks + num_predicted_peaks;
  ++num_predicted_peaks;
  return( next_predicted_peak );
}

static double twopi = 6.2831853071795865;

/*
** initializes peak values;
** some of the structure elements are no longer necessary
*/
#ifdef ANSI_C
int init_peak( Peak *peak, FLOAT pred_location, FLOAT pred_period,
               FLOAT total_signal, FLOAT proportion_fitted )
#else
int init_peak( peak, pred_location, pred_period,
               total_signal, proportion_fitted )
Peak *peak;
FLOAT pred_location, pred_period;
FLOAT total_signal, proportion_fitted;
#endif
{
  peak->pred_location = pred_location;
  peak->pred_period = pred_period;
  peak->total_signal = total_signal;
  peak->proportion_fitted = proportion_fitted;
  peak->obs_peak = 0;
  peak->best_obs_peak = 0;
  peak->best_uncalled_peak = 0;
  peak->fixed = 0;

  /*
  ** following not really necessary
  */
/*
  peak->best_nuc[0] = 4;
  peak->best_nuc[1] = 4;
  peak->best_nuc[2] = 4;
*/

  /*
  ** set values for 'N'
  */
/*
  peak->signals[4] = 0;
  peak->signal_locations[4] = peak->pred_location;
  peak->signal_type[4] = 2;
  peak->num_exposed_peaks = 0;
  peak->selected_nuc = 4;
*/
  peak->next = NULL;

  return( 0 );
}



/*
** Weight vector: weight = 1 at site,
** decreases to 0 at 0 and length.
** This is used as a triangular filter
** applied to the data before the FFT
*/
#ifdef ANSI_C
FLOAT weight_vec( FLOAT *temp_vec, FLOAT *vec,
                  int length, int site )
#else
FLOAT weight_vec( temp_vec, vec, length, site )
FLOAT *temp_vec, *vec;
int length,site;
#endif
{
  int j;
  FLOAT fac;
  FLOAT t_signal;
  static int prev_length = 0;
  static int prev_site = 0;
  static FLOAT filter[2048];
  FLOAT *pfilter;

  t_signal = 0.0;

  if( length != prev_length || site != prev_site )
  {
    if( length > 2048 )
    {
      fprintf( stderr, "weight_vec: length > array dimension\n" );
      return( 0.0 );
    }

    /*
    ** build filter
    */
    fac = 1.570796 / (FLOAT)site;
    for( j = 0; j < site; ++j )
    {
      filter[j] = (FLOAT)sin( (double)( (FLOAT)j * fac ) );
      t_signal += temp_vec[j] = filter[j] * vec[j];
    }
    fac = 1.570796 / (FLOAT)( length - site - 1 );
    for( j = site; j < length; ++j )
    {
      filter[j] = (FLOAT)cos( (double)( (FLOAT)( j - site ) * fac ) );
      t_signal += temp_vec[j] = filter[j] * vec[j];
    }
  }
  else
  {
    pfilter = filter;
    for( j = 0; j < length; ++j )
    {
      t_signal += *temp_vec = *pfilter * *vec;
      ++temp_vec;
      ++vec;
      ++pfilter;
    }
  }
  prev_length = length;
  prev_site = site;

  return( t_signal );
}



/*
**
*/
#ifdef ANSI_C
FLOAT nearest_peak( FLOAT point, FLOAT offset, FLOAT period )
#else
FLOAT nearest_peak( point, offset, period )
FLOAT point, offset, period;
#endif
{
  int i;
  FLOAT x;

  x = ( point - offset ) / period + .5;
  i = (x > 0) ? x : x - 1;
  return( offset + i * period );
}



/*
** calculate inner product of
** vec with the sin/cos of period
*/
#ifdef ANSI_C
double inner( FLOAT *vec, int length, double period, double *sum_sin )
#else
double inner( vec, length, period, sum_sin )
FLOAT *vec;
int length;
double period;
double *sum_sin;
#endif
{
  register int i;
  register double sum_cos;
  register double tmp_sin;
  register double *pcos;
  register double *psin;
  register FLOAT *pvec;
  Cos_sin_array *cs;

  cs = make_cos_sin( period, length );
  sum_cos = 0.0;
  tmp_sin = 0.0;
  pcos = cs->cos;
  psin = cs->sin;
  pvec = vec;
  for( i = 0; i < length; i++ )
  {
    sum_cos += *pvec * *pcos;
    tmp_sin += *pvec * *psin;
    ++pcos;
    ++psin;
    ++pvec;
  }    
  *sum_sin = tmp_sin;

  return( sum_cos );
}



/*
** Keeps list of cosine and sine values in arrays
** of a given length and period. When called, looks
** to see if there is an appropriate array already
** in the list; if there is, it is returned,
** otherwise a new one is created.
*/
#ifdef ANSI_C
Cos_sin_array *make_cos_sin( double period, int length )
#else
Cos_sin_array *make_cos_sin( period, length )
double period;
int length;
#endif
{
  register int i;
  register double incr_cos;
  register double incr_sin;
  register double *pcos;
  register double *psin;
  Cos_sin_array *cs;

  /* always round to nearest .001 -- to avoid mismatches
  ** due to rounding
  */
  period = ( (int)( period * 1000.0 + .5 ) ) / 1000.0;

  for( cs = alloc_cos_sin; cs; cs = cs->next )
  {
    if( period == cs->period && length <= cs->length )
    {
      return( cs );
    }
  }

  cs = (Cos_sin_array *)ourMalloc( sizeof( Cos_sin_array ) );
  cs->cos = (double *)ourMalloc( length * sizeof( double ) );
  cs->sin = (double *)ourMalloc( length * sizeof( double ) );
  cs->period = period;
  cs->length = length;
  incr_cos = cos( twopi / period );
  incr_sin = sin( twopi / period );
  pcos = cs->cos;
  psin = cs->sin;
  pcos[0] = 1.0;   
  psin[0] = 0.0;   
  for( i = 1; i < length; i++ )
  {
    pcos[i] = pcos[i-1] * incr_cos - psin[i-1] * incr_sin;
    psin[i] = psin[i-1] * incr_cos + pcos[i-1] * incr_sin;
  }    
  cs->next = alloc_cos_sin;
  alloc_cos_sin = cs;
  return( cs );
}



/*
** find_best:
**
** Purpose
** =======
**
** Find best period and the phase offset.
**
**
** Arguments
** =========
**
**         tot_vals       (input)  FLOAT array, dimension [length]
**                        Sum of traces at each point.
**
**         length         (input)  int
**                        Number of points in trace.
**
**         start_period   (input)  FLOAT
**                        First value of period to test --- a mid-value.
**
**         min_period     (input)  FLOAT
**                        Minimum value of period to test.
**
**         max_period     (input)  FLOAT
**                        Maximum value of period to test.
**
**         period_incr    (input)  FLOAT
**                        Increment (decrement) of test period
**                        from start_period
**
**         add_max_value  (output) FLOAT (*)
**                        Inner product magnitude at add_best_period.
**
**         add_best_period (output) FLOAT (*)
**                        Best period.
**
**         add_best_offset (output) FLOAT (*)
**                        Phase offset at add_best_period.
**
**         find_best      (return) int
**                        Set to 1 when stepping back provides a
**                        significantly better result.
**
*/
#ifdef ANSI_C
int find_best( FLOAT *tot_vals, int length, FLOAT start_period,
               FLOAT min_period, FLOAT max_period, FLOAT period_incr,
               FLOAT *add_max_value,
               FLOAT *add_best_period,
               FLOAT *add_best_offset )
#else
int find_best( tot_vals, length, start_period, min_period, max_period,
               period_incr, add_max_value, add_best_period, add_best_offset )
FLOAT *tot_vals;
int length;
FLOAT start_period, min_period, max_period, period_incr;
FLOAT *add_max_value, *add_best_period, *add_best_offset;
#endif
{
  int flag;
  FLOAT max_value, best_period, best_offset;
  double period;
  double temp, temp1, temp2;

  if( min_period < 2.0 )
  {
    min_period = 2.0;
  }

  max_value = 0.0;
  best_period = 0.0;
  best_offset = 0.0;

  /*
  ** Check quality of sinusoid fit starting with
  ** period = start_period and incrementing by
  ** period_incr until the quality dips, or max_period
  ** is reached.  If it never improved, look from
  ** start_period down to min_period incrementing by
  ** period_incr.
  */
  period_incr *= -1.0;
  for( period = start_period;
       period <= max_period && period >= min_period;
       period += period_incr )
  {
    /*
    ** Calculate the inner product of
    ** the tot_vals vector with
    ** the sinusoids of period.
    */
    temp1 = inner( tot_vals, length, period, &temp2 );

    /*
    ** magnitude of inner product
    */
    temp = sqrt( temp1 * temp1 + temp2 * temp2 );

    /*
    ** is this an improvement?
    */
    if( temp > max_value )
    {
      /*
      ** better, so keep the result and
      ** calculate the phase offset
      */
      max_value = temp;
      best_period = period;
      best_offset = period * acos( temp1 / temp ) / twopi;
      if( temp2 < 0.0 )
      {
        best_offset *= -1.0;
      }
    }
    else
    if( period_incr < 0.0 && best_period == start_period )
    {
      /*
      ** no improvement after first step up
      ** so start looking down from start_period.
      */
      period = start_period;
      period_incr *= -1.0;
    }
    else
    {
      /*
      ** 1. found no improvement after stepping up at least once,
      **
      ** or
      **
      ** 2. found no improvement stepping up first time so tried
      **    stepping down and found no improvement after one or more
      **    steps down.
      */
      break;
    }
  }

  *add_max_value = max_value;
  *add_best_period = best_period;
  *add_best_offset = best_offset;
  flag = 0;

  return( flag );

}



/*
** free predicted peak memory
*/
#ifdef ANSI_C
int free_predicted_peak( void )
#else
int free_predicted_peak()
#endif
{ 
  ourFree( (char *)alloc_predicted_peaks );
  return( 0 );
}   


/*
** free Cos_sin_array
*/
#ifdef ANSI_C
int free_cos_sin( void )
#else
int free_cos_sin()
#endif
{
  Cos_sin_array *cs, *next_cs;
  for( cs = alloc_cos_sin; cs; cs = next_cs )
  {
    next_cs = cs->next;
    free( cs->sin );
    free( cs->cos );
    free( cs );
  }
  alloc_cos_sin = NULL;
  return( 0 );
}
