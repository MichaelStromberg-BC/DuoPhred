/** fit_peaks.c **/

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
**    * find_observed_peaks() benefits from ideas in code written by *       **
**    * M.C. Wendl and W. Gish.                                      *       ** 
**                                                                           **
*******************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "phred.h"


/*
** make_observed_peak:
**
** Purpose
** =======
** allocates a new observed peak
** (if necessary) and
** initializes some variables
**
** Arguments
** =========
**
**         nuc:                 (input)  int
**                              Peak nucleotide.
**
**                              Nucleotide Type
**
**                              nuc     nucleotide
**                              ---     ----------
**                               0       A
**                               1       C
**                               2       G
**                               3       T
**                               4       N
**
**         make_observed_peak:  (return) Observed_peak (*)
**                              New Observed_peak structure element
**
*/
#ifdef ANSI_C
Observed_peak *make_observed_peak( int nuc )
#else
Observed_peak *make_observed_peak( nuc )
int nuc;
#endif
{
  int m, n;
  int numByte;
  Observed_peak *obs_peak;

  numByte = sizeof( Observed_peak );
  obs_peak = (Observed_peak *)ourMalloc( numByte ); 

  /*
  ** initialize structure elements
  */
  for( m = 1; m <= MAX_NUM_SPLITS; m++ )
  {
    obs_peak->peak[m] = 0;
    for( n = 1; n <= m; n++ )
    {
      obs_peak->split[m][n] = 0;
    }
  }
  for( m = 0; m <= MAX_LEFT_SHIFT + MAX_RIGHT_SHIFT; m++ )
  {
    obs_peak->dp[m] = 0;
    obs_peak->scores[m] = 0.0;
  }
  obs_peak->nuc = nuc;
  obs_peak->first_peak_index = 0;
  obs_peak->last_peak_index = 0;
  obs_peak->shift = 0;
  obs_peak->type = 0;
  obs_peak->next = 0;
  obs_peak->prev = 0;

  return( obs_peak );

} /* make_observed_peak */


/*
** insert_observed_peak:
**
** Purpose
** =======
**
** Inserts an observed peak into the doubly linked list.
** The dummy initial peak will always stay at the beginning;
** but other peaks will be sorted in order of their locations
** (defined as area median of peak)
**
** Arguments
** =========
**
**         obs_peak:              (input)  Observed_peak (*)
**                                Peak to insert into list.
**
**         insert_observed_peak:  (return) int
**                                Not meaningful.
**
*/
#ifdef ANSI_C
static int insert_observed_peak( Observed_peak *obs_peak )
#else
static int insert_observed_peak( obs_peak )
Observed_peak *obs_peak;
#endif
{
  static Observed_peak *current_peak;

  if( !obs_peak )
  {
    current_peak = 0;
  }
  else
  if( !current_peak )
  {
    /*
    ** first peak
    */
    current_peak = obs_peak;
  }
  else
  {
    /*
    ** locate insertion point
    */
    while( current_peak->next &&
           obs_peak->split[1][1] > current_peak->split[1][1] )
    {
      current_peak = current_peak->next;
    }
    while( current_peak->prev &&
           obs_peak->split[1][1] < current_peak->split[1][1] )
    {
      current_peak = current_peak->prev;
    }

    /*
    ** link structure into place
    */
    obs_peak->prev             = current_peak;
    obs_peak->next             = current_peak->next;
    if( current_peak->next )
    {
      current_peak->next->prev = obs_peak;
    }
    current_peak->next         = obs_peak;
    current_peak               = obs_peak;
  }

  return( 0 );

} /* insert_observed_peak */



/*
** find_observed_peaks:
**
** Purpose
** =======
**
** Find observed peaks.
**
** Arguments
** =========
**
**         tr_vals:              (input) FLOAT array, dimensions [nuc][point]
**                               Trace values for nucleotides 'nuc' and
**                               point 'point'.
**
**
**         tr_length:            (input) int
**                               Number of data points in trace data.
**
**
**         find_observed_peaks:  (return) Observed_peak (*)
**                               Pointer to the first element in a linked
**                               list of Observed_peaks.
**
*/
#ifdef ANSI_C
Observed_peak *find_observed_peaks( FLOAT **tr_vals, int tr_length )
#else
Observed_peak *find_observed_peaks( tr_vals, tr_length )
FLOAT **tr_vals;
int tr_length;
#endif
{
  int i, j, k, m, n, n10;
  int i_most;
  int prev_start;
  int truncFlag;
  int loPnt, hiPnt;
  int loc;
  FLOAT prev_area, peak_area, temp_area, last10, relative_area;
  FLOAT max_val[4];
  FLOAT x, last_area;
  Observed_peak *obs_peak, *first_obs_peak, *obs_peak_10_back;

  /*
  ** initialize observed peak structure
  */
  insert_observed_peak( (Observed_peak *)0 );
  first_obs_peak = make_observed_peak( -1 );
  insert_observed_peak( first_obs_peak );

  /*
  ** Get maximum trace values.
  */
  getMaxVal( max_val );

  /*
  ** locate observed peaks
  */

  /*
  ** loop through nucleotides
  */
  for( j = 0; j < 4; j++ )
  {

    /*
    ** initialize
    */
    n10 = 0;
    prev_start = 0;
    last10 = 1.0;
    prev_area = 0.0;
    truncFlag = 0;

    /*
    ** scan trace
    */
    for( i = 1; i < tr_length - 1; i++ )
    {

      /*
      ** is there a peak here?
      */
      if( tr_vals[j][i] + tr_vals[j][i] > tr_vals[j][i-1] + tr_vals[j][i+1] )
      {
        /*
        ** peak detected: initialize
        */
        peak_area = 0.0;
	prev_start = i;
        loPnt = i;

        /*
        ** sum trace values for peak area
        */
        for( ; i < tr_length - 1 &&
               tr_vals[j][i] + tr_vals[j][i] >=
               tr_vals[j][i-1] + tr_vals[j][i+1];
               i++ )
        {
          peak_area += tr_vals[j][i];
          if( tr_vals[j][i] == max_val[j] )
          {
            truncFlag = 1;
          }
        }
        hiPnt = i;

        /*
        ** find peak area relative to
        ** the mean area of the previous
        ** ten peaks
        */
        relative_area = peak_area * 10.0 / last10;

        /*
        ** is this a real peak?
        ** (ignore very small peaks)
        */
        loc = ( hiPnt + loPnt ) / 2.0;
        if( ( truncFlag && tr_vals[j][loc] > 0.05 * max_val[j] ) ||
            relative_area > MIN_PEAK_RATIO &&
            peak_area > .05 * prev_area )
        {

          if( truncFlag )
          {
            if( truncFlag > 10 )
            {
              truncFlag = 0;
            }
            else
            {
              ++truncFlag;
            }
          }

          /*
          ** keep this peak
          */
          obs_peak = make_observed_peak( j ); 

          /*
          ** initialize obs_peak_10_back
          */
          if( !n10 )
          {
            obs_peak_10_back = obs_peak;
          }

          /*
          ** sum of last 10 peak areas
          */
          last10 += peak_area;

          /*
          ** detected more than 10 peaks
          */
          n10++;
          if( n10 > 10 )
          {
            /*
            ** reduce last10 by area of 11th peak back
            */
            last10 -= obs_peak_10_back->area;

            /*
            ** link in new obs_peak_10_back
            */
            for( obs_peak_10_back = obs_peak_10_back->next;
                 obs_peak_10_back->nuc != j;
                 obs_peak_10_back = obs_peak_10_back->next );
          }

          /*
          ** split each peak area MAX_NUM_SPLITS different ways:
          **
          ** split number (m)       divisions (n)
          ** ================       ===========================
          **      1                  0.5
          **      2                  0.25, 0.75
          **      3                  0.167, 0.5, 0.833
          **      4                  0.125, 0.375, 0.625, 0.875
          **      ...                ...
          **    MAX_NUM_SPLITS       ...
          **
          ** Note: There may be more divisions than data points
          **       that define the peak.
          */
          temp_area = 0.0;
          i_most = prev_start;
          for( k = prev_start; k < i; k++ )
          {
            last_area = temp_area;
            temp_area += tr_vals[j][k];
            if( tr_vals[j][k] >= tr_vals[j][i_most] )
            {
              i_most = k;
            }
            for( m = 1; m <= MAX_NUM_SPLITS; m++ )
            {
              for( n = 1; n <= m; n++ )
              {
                x = ( n + n - 1 ) * peak_area / ( m + m );
                if( temp_area < x )
                {
                  break;
                }
                else
                if( !obs_peak->split[m][n] )
                {
                  /*
                  ** store the split location ---
                  ** the data point where the division
                  ** occurs
                  */
                  obs_peak->split[m][n] = k - ( temp_area - x > x - last_area );
                }
              }
            }
          } /* for k */

          /*
          ** fill in peak structure elements and
          ** insert it into the linked list
          */
          obs_peak->i_left = prev_start;
          obs_peak->i_rite = i;
          obs_peak->i_maxx = i_most;
          obs_peak->area = peak_area;
          prev_area = peak_area;
          obs_peak->relative_nuc_area = relative_area;
          insert_observed_peak( obs_peak );

        } /* if relative_area */

      } /* if tr_vals */

    } /* for i */

  } /* for j */


  return( first_obs_peak );

} /* find_observed_peaks */
    
    

/*
** free observed peak memory
*/ 
#ifdef ANSI_C
int free_obs_peak( Observed_peak *first_obs_peak )
#else
int free_obs_peak( first_obs_peak )
Observed_peak *first_obs_peak;
#endif
{
  Observed_peak *obs_peak, *next_obs_peak;

  obs_peak = first_obs_peak;
  while( obs_peak )
  {
    next_obs_peak = obs_peak->next;
    free( obs_peak );
    obs_peak = next_obs_peak;
  }

  return( 0 );
}


/*
** fit_peaks:
**
** Purpose
** =======
**
** Call bases.
**
** Arguments
** =========
**
**         tr_length:      (input) int
**                         Number of points in trace.
**
**         tr_vals:        (input) FLOAT
**                         Traces.
**
**         first_peak:     (input) Peak (*)
**                         Linked list of predicted peaks.
**
**         first_obs_peak: (output) Observed_peak (*)
**                         Linked list of observed peaks.
**
**         fit_peaks:      (return) int
**                         Number of peaks 'fixed'.
**
*/
#ifdef ANSI_C
int fit_peaks( int tr_length, FLOAT **tr_vals, Peak *first_peak,
               Observed_peak *first_obs_peak )
#else
int fit_peaks( tr_length, tr_vals, first_peak, first_obs_peak )
int tr_length;
FLOAT **tr_vals;
Peak *first_peak;
Observed_peak *first_obs_peak;
#endif
{
  int i, j, k, m;
  int i_peak, n_peaks, max_m;
  int location, counter, n10, n_fixed;
  int first_peak_index, last_peak_index;
  int max_label;
  FLOAT score, max_score, last10;
  FLOAT shift, next_shift;
  FLOAT *score_vec_1, *score_vec_2;
  FLOAT *prev_scores, *current_scores, *temp_scores;
  int *peak_indices;
  Peak *peak, *last_peak, *best_neighbor, *peak_10_back;
  Peak **peaks;
  Observed_peak *obs_peak, *last_obs_peak, *obs_peak_10_back;

  /*
  ** Count predicted peaks and
  ** identify last predicted peak.
  */
  for( peak = first_peak, n_peaks = 0; peak; peak = peak->next )
  {
    n_peaks++;
    last_peak = peak;
  }
    
  /*
  ** Allocate array to store predicted peak pointers and scoring arrays.
  */
  peaks = (Peak **)ourMalloc( n_peaks * sizeof( Peak * ) );
  score_vec_1 = (FLOAT *)ourMalloc( n_peaks * sizeof( FLOAT ) );
  score_vec_2 = (FLOAT *)ourMalloc( n_peaks * sizeof( FLOAT ) );
 
  /*
  ** Store predicted peak pointers in array and
  ** initialize scoring vectors.
  */
  for( peak = first_peak, i_peak = 0; peak; peak = peak->next, i_peak++ )
  {
    peaks[i_peak] = peak;
    score_vec_1[i_peak] = 0.0;
    score_vec_2[i_peak] = 0.0;
  }
  
  /*
  ** Allocate peak_indices array and initialize it
  ** such that each element corresponds to a trace
  ** point and the element value is the index of the
  ** predicted peak immediately to the left.
  */
  peak_indices = (int *)ourMalloc( ( tr_length + 10 ) * sizeof( int ) );
  peak_indices[0] = 0;
  i_peak = 0;
  for( i = 1; i < tr_length + 10; i++ )
  {
    if( i_peak < n_peaks && i >= peaks[i_peak]->pred_location )
    {
      peak_indices[i] = i_peak;
      i_peak++;
    }
    else
    {
      peak_indices[i] = peak_indices[i-1];
    }
  } /* for i */
  
  /*
  ** Preliminary match of observed peaks
  ** to predicted peaks working from front to back
  ** Pointers to observed peaks are stored in
  ** best_obs_peak element of predicted peak structures.
  ** Some predicted peaks may be unmatched.
  */
  for( obs_peak = first_obs_peak->next; obs_peak; obs_peak = obs_peak->next )
  {

    /*
    ** location => location on trace of current observed peak division
    ** peak     => pointer to predicted peak at point [location]
    **             on the trace
    ** shift    => amount observed peak is displaced from predicted
    **             location, expressed as a fraction of the predicted
    **             period (a negative shift value indicates that the
    **             location lies to the right of the predicted peak).
    */
    location = obs_peak->split[1][1];
    peak = peaks[peak_indices[location]];
    shift = measure_shift( (FLOAT)location, peak );

    /*
    ** Does obs_peak lie too far to the right of
    ** the predicted peak?
    */
    if( shift <= -0.6 && peak->next )
    {
      /*
      ** Large shift, skip
      ** to next predicted peak.
      */
      peak = peak->next;
    }

    /*
    ** Store observed peak pointer in
    ** best_obs_peak element of predicted
    ** peak, if there's no best_obs_peak assigned, or
    ** if the area of the obs_peak is greater
    ** than the area of the previously assigned
    ** best_obs_peak.  Thus we find the largest obs_peak that is
    ** located in the region from about halfway between the predicted
    ** peak and the previous predicted peak and the predicted peak
    ** and the next predicted peak.
    */
    if( !peak->best_obs_peak || obs_peak->area > peak->best_obs_peak->area )
    {
      peak->best_obs_peak = obs_peak;
    }

    /*
    ** pointer to last peak in
    ** obs_peak list, used later
    */
    last_obs_peak = obs_peak;

  }  /* for obs_peak */
  

  /*
  ** Working from front to back, express area of each
  ** peak as a fraction of the mean area of the ten
  ** preceding best_obs_peaks and store in relative_area
  ** element of observed peak structure.
  */
  last10 = 1.0;
  n10 = 0;
  for( obs_peak = first_obs_peak->next; obs_peak; obs_peak = obs_peak->next )
  {
    obs_peak->relative_area = 10.0 * obs_peak->area / last10;
    peak = peaks[peak_indices[obs_peak->split[1][1]]];
    if( obs_peak == peak->best_obs_peak || peak->next &&
        obs_peak == peak->next->best_obs_peak )
    {
      if( !n10 )
      {
        obs_peak_10_back = obs_peak;
      }
      last10 += obs_peak->area;
      n10++;
      if( n10 > 10 )
      {
        last10 -= obs_peak_10_back->area;
        for( obs_peak_10_back = obs_peak_10_back->next;
             obs_peak_10_back;
             obs_peak_10_back = obs_peak_10_back->next )
        {
          peak = peaks[peak_indices[obs_peak_10_back->split[1][1]]];
          if( obs_peak_10_back == peak->best_obs_peak || peak->next &&
              obs_peak_10_back == peak->next->best_obs_peak )
          {
            break;
          }
        }
      } /* if n10 */
    } /* if obs_peak */
  } /* for obs_peak */
  
  /*
  ** Set next_no_observed flag of predicted peak to zero
  ** if either there is no observed peak assigned to the
  ** predicted peak or the relative area of the assigned
  ** observed peak is less than a minimum value.
  */
  for( peak = first_peak; peak; peak = peak->next )
  {
    peak->next_no_observed = peak->best_obs_peak &&
                             peak->best_obs_peak->relative_area >
                             MIN_RELATIVE_AREA;
  }
  

  /*
  ** Working from back to front, fill next_no_observed element
  ** of each predicted peak with an integer count of the
  ** number of consecutive valid peaks following it in the
  ** predicted peak list.  The value of the next_no_observed
  ** element of the "invalid" predicted peak is zero.
  */
  counter = 0;
  for( peak = last_peak; peak; peak = peak->prev )
  {

    /*
    ** Increment counter if
    **
    **    (1) a best_obs_peak is assigned to predicted peak i
    **    (2) a best_obs_peak is assigned to predicted peak i+1
    **    (3) the distance between the best_obs_peak assigned to
    **        predicted peak(i) and the best_obs_peak assigned to
    **        predicted peak(i+1) is within 20% of the predicted
    **        period.
    **
    ** otherwise reset counter to zero.
    */
    counter = peak->next_no_observed &&
              peak->next &&
              peak->next->next_no_observed &&
              fabs( peak->next->best_obs_peak->split[1][1] -
                    ( peak->best_obs_peak->split[1][1] +
                      peak->pred_period ) ) / peak->pred_period < 0.2 ?
                                              counter + 1 : 0;

    /*
    ** Good peak(s)?
    */
    if( counter >= 1 )
    {
      /*
      ** Set next_no_observed element
      ** of predicted peak to number
      ** of consecutive "good" peaks.
      */
      peak->next_no_observed = counter;
    }

    /*
    ** Is this the first consecutive peak (and
    ** not the end of the predicted peak list?)
    */
    if( counter == 1 && peak->next )
    {
      /*
      ** Yes, set next_no_observed element of
      ** next predicted peak to zero. (It
      ** may have been set to one in the
      ** previous loop).
      */
      peak->next->next_no_observed = 0;
    }

  } /* for peak */
  

  /*
  ** Working from front to back, assign 'fixed'
  ** observed peaks to predicted peaks.  These
  ** are 'easy calls' in that the observed peaks
  ** align well with the predicted peaks: they're
  ** considered unambiguous by phred so these
  ** calls are not reconsidered at any later time.
  */
  n_fixed = 0;
  for( peak = first_peak; peak; )
  {

    /*
    ** Are 4 or more subsequent predicted
    ** peaks matched with observed peaks?
    */
    if( peak->next_no_observed >= 4 )
    {

      /*
      ** Yes, walk down band of consecutive matched
      ** predicted peaks marking the predicted peaks
      ** as satisfied with the best_obs_peak. (No
      ** nearby compressions to make these matches
      ** questionable.)
      */
      for( ; peak && peak->next_no_observed; peak = peak->next )
      {
        peak->obs_peak = peak->best_obs_peak;
        peak->obs_peak->peak[1] = peak;
        peak->obs_peak->shift = measure_shift(
                                  (FLOAT)peak->obs_peak->split[1][1], peak );
        peak->obs_peak->type = 1;
        peak->fixed = 1;
        n_fixed++;
      }
    }
    else
    {
      peak = peak->next;
    }

  } /* for peak */


  /*
  ** Use Smith-Waterman-based comparison
  ** algorithm to consider for each position
  **   (i)  "insertion" (actually, deletion) of observed peak; i.e.
  **        the observed peak is not assigned to a predicted peak,
  **   (ii) alignment of observed peak against 1, 2, ..., MAX_NUM_SPLITS
  **        predicted peaks. For each of these, compute score.
  **
  ** Start at second observed peak and work toward end of obs_peak list.
  */
  prev_scores    = score_vec_1;
  current_scores = score_vec_2;
  
  for( obs_peak = first_obs_peak->next; obs_peak; obs_peak = obs_peak->next )
  {

    /*
    ** Begin with some setup and testing of
    ** peak 'quality'...
    **
    ** Get index to predicted peak MAX_RIGHT_SHIFT
    ** predicted peaks to left of the predicted
    ** peak nearest the current observed peak center.
    */
    first_peak_index = peak_indices[obs_peak->split[1][1]] - MAX_RIGHT_SHIFT;
    if( first_peak_index < 0 )
    {
      first_peak_index = 0;
    }
    obs_peak->first_peak_index = first_peak_index;

    /*
    ** Get index to predicted peak MAX_LEFT_SHIFT
    ** predicted peaks to right of the predicted
    ** peak nearest the current observed peak center.
    */
    last_peak_index = peak_indices[obs_peak->split[1][1]] + MAX_LEFT_SHIFT;
    if( last_peak_index >= n_peaks )
    {
      last_peak_index = n_peaks - 1;
    }
    obs_peak->last_peak_index = last_peak_index;

    /*
    ** Is this a peak from an "uncompressed" band.
    */
    if( obs_peak->type )
    {
      /*
      ** Yes, is already assigned to a peak: reset
      **  alignment score vectors within current window.
      */
      for( i_peak = first_peak_index; i_peak <= last_peak_index; i_peak++ )
      {
        current_scores[i_peak] = 0.0;
        prev_scores[i_peak]    = 0.0;
      }

      /*
      ** Skip to next obs_peak.
      */
      continue; 

    } /* if obs_peak->type */

    /*
    ** Is the relative area of this observed peak small?
    */
    if( obs_peak->relative_area < 0.1 )
    {
      /*
      ** Yes, skip to next obs_peak leaving
      ** the score vectors untouched.
      */
      continue;
    }

    /*
    ** Reset some alignment score vector elements.
    */
    for( i_peak = obs_peak->prev->first_peak_index;
         i_peak < first_peak_index;       i_peak++ )
    {
      /*
      ** Set to zero scores from the beginning of the
      ** window for the previous peak to the beginning
      ** of the window for the current peak.
      */
      current_scores[i_peak] = 0.0;
    }

    /*
    ** OK, we are now set up to begin testing for insertions
    ** and alignments...so we loop through predicted peaks
    ** within the window around the observed peak in order
    ** to calculate the scores for these various alignment
    ** possibilities (the alignments that produce these scores
    ** are encoded in the obs_peak->dp[] array).
    **
    ** i_peak corresponds to the column index of
    ** the Smith-Waterman matrix and obs_peak
    ** corresponds to the row index of the matrix
    **
    ** As a start, think this through in the case of
    ** no splitting, i.e. max_m = 1, m = 1, and k = 1.
    */
    for( i_peak = first_peak_index; i_peak <= last_peak_index; i_peak++ )
    {

      /*
      ** this score corresponds to deletion
      ** of current (observed) peak; it is
      ** the score for this predicted peak
      ** as computed for the previous observed
      ** peak (row of the Smith-Waterman matrix)
      */
      max_score = prev_scores[i_peak];
      max_label = 0;

      /*
      ** let's see what happens when
      ** we align the observed peak
      ** split in 1, 2, ..., MAX_NUM_SPLITS
      ** divisions with an equal number of
      ** of predicted peaks.
      */
      if( obs_peak->relative_area < MIN_SPLITTABLE_AREA )
      {
        /*
        ** don't split peaks that don't exceed minimum
        ** relative peak area
        */
        max_m = 1;
      }
      else
      {
        /*
        ** split peak into MAX_NUM_SPLITS
        ** (or fewer at start of window)
        */
        max_m = ( i_peak + 1 > MAX_NUM_SPLITS ) ? MAX_NUM_SPLITS : i_peak + 1;
      }
      
      /*
      ** look at each way of dividing
      ** the peak area (the number of
      ** divisions)
      */
      for( m = 1; m <= max_m; m++ )
      {

        /*
        ** initial score is the previous score "m"
        ** predicted peaks to the left of the current
        ** predicted peak.  This is an extension of the
        ** alignment; that is, the
        ** H(i,j) = H(i-1,j-1) + score(i,j)
        ** of Smith-Waterman.
        */
        score = i_peak >= m ? prev_scores[i_peak - m] : 0.0;

        /*
        ** do the divisions of the split peak
        ** align with the set of corresponding
        ** predicted peaks?
        */
        for( k = 1; k <= m; k++ )
        {

          /*
          ** the corresponding predicted peak
          ** for this position in the split
          */
          peak = peaks[i_peak - m + k];

          /*
          ** measure the difference between the positions
          ** of the predicted peak and the "area division" of
          ** the current observed peak.  Express the difference
          ** as a fraction of the distance between predicted peaks
          ** (the period).
          */
          shift = measure_shift( (FLOAT)obs_peak->split[m][k], peak );

          /*
          ** ignore this split division if the shift is too large
          */
          if( shift < -MAX_RIGHT_DISPLACEMENT || shift > MAX_LEFT_DISPLACEMENT )
          {
            break;
          }
          score += alignment_score( obs_peak, m, k, peak, shift );

        } /* for k */

        /*
        ** record a maximum score -- if we didn't
        ** break out of the "k" loop prematurely;
        ** i.e., the split division wasn't
        ** clearly bad.
        */
        if( k > m && score > max_score )
        {
          max_score = score;
          max_label = m;
        }

      } /* for m */

      /*
      ** record the results of our alignment
      ** tests
      */
      current_scores[i_peak] = max_score;
      obs_peak->dp[i_peak - first_peak_index] = max_label;
      obs_peak->scores[i_peak - first_peak_index] = max_score;

    } /* for i_peak */

    /*
    ** exchange current and previous scores
    */
    temp_scores = prev_scores;
    prev_scores = current_scores;
    current_scores = temp_scores;

  } /* for obs_peak */

  /*
  ** last peak unassigned -- but should be filled in later
  */
  i_peak = n_peaks - 2;


  /*
  ** OK...we have the scores for the best alignments,
  ** now we reconstruct good observed/predicted peak
  ** alignment by backtracking...
  */

  /* We loop through the observed peaks beginning with
  ** the last observed peak and working toward the beginning
  ** in order to assign the observed peaks to predicted
  ** peaks.
  **
  ** The essential information used in backtracking is
  **
  ** a) the dynamic programming vector assigned to each
  **    observed peak (obs_peak->dp) and the dynamic
  **    programming score vector assigned to each observed
  **    peak (obs_peak->scores).  The dp[] vector describes
  **    whether/how to split the peak and the scores[]
  **    vector indicates whether or not the observed peak
  **    may be assigned to the predicted peak corresponding
  **    to the vector index.  The obs_peak->first_peak_index
  **    correlates the vector elements to the predicted peak
  **    indices.
  **
  ** b) whether an observed peak was assigned to a predicted
  **    peak already.  Observed peak were assigned to predicted
  **    peaks on creation of "fixed" peaks.  This was the only
  **    assignment of observed peaks to predicted peaks before
  **    this next loop,
  **
  ** c) the "type" element of each observed peak, which
  **    was initialized to zero on creation of the
  **    observed peak.  At this point only "fixed" peaks
  **    have non-zero values of the element "type", and
  **    those are set to one,
  **
  ** d) the "fixed" element of each predicted peak, which
  **    was initialized to zero on creation of the
  **    observed peak.  The value of "fixed" for each
  **    fixed peak was set to one.  These values are
  **    untouched from this point onward.
  **
  */

  for( obs_peak = last_obs_peak;
       obs_peak && i_peak >= 0;
       obs_peak = obs_peak->prev )
  {

    /*
    ** Work from the current observed peak toward the
    ** front of the list until we find an unmatched
    ** observed peak, and setting i_peak to the nearest
    ** predicted peak, if that observed peak is not
    ** already assigned to the predicted peak; otherwise
    ** select the previous predicted peak.  So the selected
    ** predicted peak is within one predicted period and
    ** to the left of the observed peak or it is within
    ** one predicted period and to the right of the
    ** of the observed peak.
    */
    for( ; obs_peak && obs_peak->type; obs_peak = obs_peak->prev )
    {
      i = peak_indices[obs_peak->split[1][1]];
      i_peak = ( peaks[i]->obs_peak == obs_peak ) ? i - 1 : i;
    }

    /*
    ** Step toward the beginning of the predicted peaks.  On
    ** finding a "fixed" predicted peak, reset the observed peak
    ** to the one immediately to the left of the observed
    ** peak assigned to the "fixed" predicted peak, and continue.
    ** Stop when (1) the predicted peak location is to the left of
    ** the right edge of the shift window (defined between
    ** first_peak_index and last_peak_index) and (2) the difference
    ** between the shift for this observed/predicted peak pair does
    ** not exceed the shift of the obs_peak assigned to the next
    ** predicted peak by 0.7.
    */
    for( ; obs_peak && i_peak >= 0; i_peak-- )
    {
      /*
      ** Is this predicted peak matched (fixed)?
      */
      if( peaks[i_peak]->fixed )
      {
        /*
        ** Yes, move left to next observed peak.
        */
        obs_peak = peaks[i_peak]->obs_peak->prev;
      }
      else
      if( i_peak <= obs_peak->last_peak_index )
      {
        /*
        ** The predicted peak is not outside
        ** the right edge of the window for the
        ** observed peak.
        */
        shift = measure_shift( (FLOAT)obs_peak->split[1][1], peaks[i_peak] );
	next_shift = peaks[i_peak+1]->obs_peak ?
                     peaks[i_peak+1]->obs_peak->shift : 0.0;

        /*
        ** Accept this predicted peak, shift change
        ** is within acceptance range.  This breaks
        ** out of the predicted peak stepping loop
        ** when the selected predicted peak becomes
        ** unreasonably distant from the observed peak,
        ** in which case, we try the next obs_peak and
        ** reset the selected predicted peak back down
        ** next to the new observed peak.
        */
        if( shift - next_shift < MAX_SHIFT_CHANGE )
        {
          break;
        }

      } /* if i_peak */

    } /* for ; obs_peak */

    if( !obs_peak || i_peak < 0 )
    {
      break;
    }

    if( obs_peak->type )
    {
      continue;
    }

    if( shift < -0.4 &&
        peaks[i_peak+1]->obs_peak && 
        peaks[i_peak+1]->obs_peak->shift > 0.4 )
    {
      /*
      ** peaks that are too close together
      ** (and not in compression); need to make this
      ** condition more flexible
      */
      continue;
    } /* if shift */

    /*
    ** Assign the observed peak to the predicted peak
    ** if the dp score is non-zero.  (Apply the determined
    ** splitting when m > 1.)
    */
    if( obs_peak->scores[i_peak - obs_peak->first_peak_index] &&
	(m  = obs_peak->dp[i_peak - obs_peak->first_peak_index] ) )
    {
      obs_peak->type = m;
      if( peaks[i_peak]->obs_peak )
      {
        peaks[i_peak]->obs_peak->type = 0;
      }

      /*
      ** Match the (split) observed peak and
      ** corresponding predicted peak(s).
      */
      for( j = m; j > 0 && i_peak >= 0; j-- )
      { 
	peaks[i_peak]->obs_peak = obs_peak;
	obs_peak->peak[j] = peaks[i_peak];
	obs_peak->shift = measure_shift( (FLOAT)obs_peak->split[m][j],
                                         peaks[i_peak] );
	i_peak--;
      }
    } /* if obs_peak->scores[] */
  } /* for obs_peak */

  /*
  ** Loop through observed peaks in order to assign
  ** observed peaks to best_obs_peak and best_uncalled_peak
  ** elements of predicted peaks.
  */
  last10 = 1.0;
  n10 = 0;
  for( obs_peak = first_obs_peak->next; obs_peak; obs_peak = obs_peak->next )
  {

    /*
    ** Record corrected relative areas of observed peaks including
    ** those split in the previous loop.
    */
    obs_peak->relative_area = 10.0 * obs_peak->area / last10;
    if( obs_peak->type > 1)
    {
      obs_peak->relative_area /= obs_peak->type;
    }

    /*
    ** Is observed peak matched to predicted peak?
    */
    if( obs_peak->type )
    {
      /*
      ** Yes - update the last ten peak area sum.  Then
      ** continue with next obs_peak.
      */
      if( !n10 )
      {
        obs_peak_10_back = obs_peak;
      }
      last10 += obs_peak->area / obs_peak->type;
      n10++;
      if( n10 > 10 )
      {
        last10 -= obs_peak_10_back->area / obs_peak_10_back->type;
        for( obs_peak_10_back = obs_peak_10_back->next;
             !obs_peak_10_back->type;
             obs_peak_10_back = obs_peak_10_back->next );
      }
    }
    else
    {
      /*
      ** No - try a tentative match of predicted and observed peak.
      */
      location = obs_peak->split[1][1]; 
      peak = peaks[peak_indices[location]];

      if( measure_shift( (FLOAT)location, peak ) <= -0.6 && peak->next )
      {
        peak = peak->next;
      }
      obs_peak->shift = measure_shift( (FLOAT)obs_peak->split[1][1], peak );

      /*
      ** Is an observed peak assigned to this predicted peak?
      */
      if( !peak->obs_peak )
      {
        /*
        ** No observed peak for this predicted peak.
        */
        if( !peak->best_obs_peak || peak->best_obs_peak->type )
        {
          /*
          ** No best_obs_peak assigned or best_obs_peak is called already,
          ** so assigned this obs_peak to best_obs_peak.
          */
          peak->best_obs_peak = obs_peak;
        }
        else
        if( obs_peak->relative_area > peak->best_obs_peak->relative_area )
        {
          if( obs_peak->nuc != peak->best_obs_peak->nuc )
          {
            peak->best_uncalled_peak = peak->best_obs_peak;
          }
          peak->best_obs_peak = obs_peak;
        }
        else
        if( !peak->best_uncalled_peak ||
            obs_peak->relative_area >
            peak->best_uncalled_peak->relative_area )
        {
          if( obs_peak->nuc != peak->best_obs_peak->nuc ||
            ( tr_vals[obs_peak->nuc][(obs_peak->split[1][1]+peak->best_obs_peak->split[1][1])/2] <
              0.95 * tr_vals[obs_peak->nuc][obs_peak->split[1][1]] &&
              tr_vals[obs_peak->nuc][(obs_peak->split[1][1]+peak->best_obs_peak->split[1][1])/2] <
              0.95 * tr_vals[obs_peak->nuc][peak->best_obs_peak->split[1][1]] ) )
          {
            peak->best_uncalled_peak = obs_peak;
          }
        }
      }
      else
      if( !peak->best_uncalled_peak ||
          obs_peak->relative_area >
          peak->best_uncalled_peak->relative_area )
      {
        if( obs_peak->nuc != peak->obs_peak->nuc ||
            ( tr_vals[obs_peak->nuc][(obs_peak->split[1][1]+peak->obs_peak->split[1][1])/2] <
              0.95 * tr_vals[obs_peak->nuc][obs_peak->split[1][1]] &&
              tr_vals[obs_peak->nuc][(obs_peak->split[1][1]+peak->obs_peak->split[1][1])/2] <
              0.95 * tr_vals[obs_peak->nuc][peak->obs_peak->split[1][1]] ) )
        {
          peak->best_uncalled_peak = obs_peak;
        }
      }
    } /* if obs_peak->type */
  } /* for obs_peak */

  /*
  ** Loop through predicted peaks trying to assign
  ** observed peaks using either a best_obs_peak or
  ** by splitting neighboring peak.
  */
  for( peak = first_peak->next; peak; peak = peak->next )
  {
    if( !peak->obs_peak )
    {
      if( peak->best_obs_peak && !peak->best_obs_peak->type )
      {
        if( peak->next &&
            peak->next->obs_peak &&
            peak->next->obs_peak->type == 1 &&
            peak->best_obs_peak->split[1][1] > peak->next->obs_peak->split[1][1] )
        {
            peak->obs_peak = peak->next->obs_peak;
            peak->obs_peak->peak[1] = peak;
            peak->obs_peak->shift = measure_shift(
                                           (FLOAT)peak->obs_peak->split[1][1],
                                           peak );
            peak->next->obs_peak = peak->best_obs_peak;
            peak->next->obs_peak->type = 1;
            peak->next->obs_peak->peak[1] = peak->next;
            peak->next->obs_peak->shift = measure_shift(
                                           (FLOAT)peak->next->obs_peak->split[1][1],
                                           peak );
        }
        else
        {
            peak->obs_peak = peak->best_obs_peak;
            peak->obs_peak->type = 1;
            peak->obs_peak->peak[1] = peak;
            peak->obs_peak->shift = measure_shift(
                                           (FLOAT)peak->obs_peak->split[1][1],
                                           peak );
        }
      }
      else
      {
        Peak *neighbor_1, *neighbor_2;
        neighbor_1 = 0;
        neighbor_2 = 0;
        if( peak->prev &&
            peak->prev->obs_peak &&
            peak->prev->obs_peak->type < MAX_NUM_SPLITS &&
            peak->prev->obs_peak->relative_area > MIN_SPLITTABLE_AREA )
        {
          neighbor_1 = peak->prev;
        }
        if( peak->next &&
            peak->next->obs_peak &&
            peak->next->obs_peak->type < MAX_NUM_SPLITS &&
            peak->next->obs_peak->relative_area > MIN_SPLITTABLE_AREA )
        {
          neighbor_2 = peak->next;
        }
        if( neighbor_1 && neighbor_2 )
        {
          if( neighbor_1->obs_peak->relative_area >
              neighbor_2->obs_peak->relative_area )
          {
            best_neighbor = neighbor_1;
          }
          else
          {
            best_neighbor = neighbor_2;
          }
        }
        else
        if( neighbor_1 )
        {
          best_neighbor = neighbor_1;
        }
        else
        if( neighbor_2 )
        {
          best_neighbor = neighbor_2;
        }
        else
        {
          best_neighbor = 0;
        }

        if( best_neighbor )
        {
          peak->obs_peak = best_neighbor->obs_peak;
          peak->obs_peak->type += 1;
          peak->obs_peak->peak[peak->obs_peak->type] = peak;
          if( peak->obs_peak->type > 1 )
          {
            peak->obs_peak->relative_area *= (FLOAT)( peak->obs_peak->type - 1 ) /
                                             (FLOAT)( peak->obs_peak->type );
          }
        }
      } /* peak->best_obs_peak */

    } /* if !peak->obs_peak */
    else
    if( peak->best_obs_peak &&
        !peak->best_obs_peak->type &&
        ( !peak->best_uncalled_peak ||
          peak->best_obs_peak->relative_area >
          peak->best_uncalled_peak->relative_area ) )
    {
      peak->best_uncalled_peak = peak->best_obs_peak;
    }

  } /* for peak */

  /*
  ** Update relative peak areas.
  */
  for( peak = first_peak; peak; peak = peak->next )
  {
    if( peak->obs_peak )
    {
      break;
    }
  }
  n10 = 1;
  obs_peak = peak->obs_peak;
  last10 = obs_peak->area / (float)obs_peak->type;
  obs_peak->relative_area = 1.0;
  peak_10_back = peak;
  for( peak = peak->next; peak; peak = peak->next )
  {
    obs_peak = peak->obs_peak;

    if( obs_peak )
    {
      obs_peak->relative_area = ( obs_peak->area / (float)obs_peak->type ) /
                                ( last10         / (float)n10 );
      last10 += obs_peak->area / (float)obs_peak->type;

      if( n10 >= 10 )
      {
        while( peak_10_back->obs_peak == NULL )
        {
          peak_10_back = peak_10_back->next;
        }
        last10 -= peak_10_back->obs_peak->area / (float)peak_10_back->obs_peak->type;
        peak_10_back = peak_10_back->next;
      }
      else
      {
        ++n10;
      }
    }
  }

  ourFree( (char *)peaks );
  ourFree( (char *)peak_indices );
  ourFree( (char *)score_vec_1 );
  ourFree( (char *)score_vec_2 );

  return( n_fixed );

} /* fit_peaks */



/*
** alignment_score:
**
** Purpose
** =======
**
** calculate alignment score.
**
** Arguments
** =========
**
**         obs_peak         (input) Observed_peak (*)
**                          Pointer to Observed_peak list.
**
**         m                (input) int
**                          Unused.
**
**         k                (input) int
**                          Number of splits for peak.
**
**         peak             (input) Peak (*)
**                          Unused.
**
**         shift            (input) FLOAT
**                          Shift of obs_peak from predicted
**                          position, expressed as a fraction
**                          of the distance between peaks.
**
**         alignment_score  (return) FLOAT
**                          Alignment score.
**
*/
#ifdef ANSI_C
FLOAT alignment_score( Observed_peak *obs_peak, int m,
                       int k, Peak *peak, FLOAT shift )
#else
FLOAT alignment_score( obs_peak, m, k, peak, shift )
Observed_peak *obs_peak;
int m, k;
Peak *peak;
FLOAT shift;
#endif
{
  int i;
  int n_offsets;
  FLOAT area, score;
  FLOAT right_shift_penalty, left_shift_penalty;

  /*
  ** note sliding area for split peaks
  */
  area = k == 1 ? obs_peak->area :
                  ( 1.0 - SPLIT_PENALTY ) * obs_peak->area / k;

  if( k == 1 )
  {
    left_shift_penalty = LEFT_SHIFT_PENALTY; 
    right_shift_penalty = RIGHT_SHIFT_PENALTY;
  }
  else
  {
    left_shift_penalty = LEFT_SPLIT_SHIFT_PENALTY; 
    right_shift_penalty = RIGHT_SPLIT_SHIFT_PENALTY;
  }

  score = area * 0.01;
  if( shift < 0.0 )
  {
    /*
    ** score for right shifted peak
    */
    score *= 1.0 + shift * right_shift_penalty;
  }
  else
  {
    /*
    ** score for left shifted peak
    */
    n_offsets = shift;
    for( i = 0; i < n_offsets; i++ )
    {
      score *= ( 1.0 - left_shift_penalty );
    }
    score *= 1.0 - ( shift - n_offsets ) * left_shift_penalty;
  }

   return( score );
  
} /* alignment_score */



/*
** measure_shift:
**
** Purpose
** =======
**
** Measure shift of observed peak from the predicted peak as a
** fraction of the distance between predicted peaks.
**
** Arguments
** =========
**
**         location    (input)  FLOAT
**                     Location of observed peak.
**
**         peak        (input)  Peak (*)
**                     Pointer to a predicted peak.
**
**         measure_shift  (return)  FLOAT
**                        Observed peak shift form predicted position
**                        expressed as a fraction of the distance
**                        between two adjacent predicted peaks.
**                        A positive value indicates that the observed
**                        peak is left of the predicted peak (a possible
**                        compression).
**
*/
#ifdef ANSI_C
FLOAT measure_shift( FLOAT location, Peak *peak )
#else
FLOAT measure_shift( location, peak )
FLOAT location;
Peak *peak;
#endif
{
  FLOAT shift, pred_loc;

  pred_loc = peak->pred_location;

  /*
  ** Is the observed peak located to the right of the
  ** predicted peak?
  */
  if( location > pred_loc )
  {
    /*
    ** Calculate the shift relative to the distance between
    ** this predicted peak and the next predicted peak to
    ** the right.
    */
    shift = peak->next ?
            ( pred_loc - location ) /
            ( peak->next->pred_location - pred_loc ) : 0;
  }
  else
  {
    /*
    ** Calculate the shift relative to the distance between
    ** this predicted peak and the next predicted peak to
    ** the left.
    */
    shift = peak->prev ?
            ( pred_loc - location ) /
            ( pred_loc - peak->prev->pred_location ) : 0;
  }

  return( shift );

}

