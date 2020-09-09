/** fit_compress.c **/

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
int fit_compress( Peak *first_peak, int tr_length, FLOAT **tr_vals, FLOAT *tot_vals )
#else
int fit_compress( first_peak, tr_length, tr_vals, tot_vals )
Peak *first_peak;
int tr_length;
FLOAT **tr_vals;
FLOAT *tot_vals;
#endif
{
  int i, j, n;
  int shift;
  int type;
  int numPeak;
  int numByte;
  int *inx;
  Observed_peak *obs_peak;
  Peak *peak;
  Peak *new_peak;
  char *seq;
  char base[4];
  FLOAT fac;
  FLOAT *relative_area;
  FLOAT *area;
  FLOAT *area_sum;
  FLOAT *area_var;
  int *compress;
  BaseQual *baseQual;

  /*
  ** Count peaks.
  */
  numPeak = 0;
  peak = first_peak;
  while( peak )
  {
    ++numPeak;
    peak = peak->next;
  }

  /*
  ** Allocate memory.
  */
  numByte = numPeak * sizeof( char );
  seq = (char *)malloc( numByte );
  if( seq == NULL )
  {
    fprintf( stderr, "fit_compress: error: unable to allocate memory\n" );
    return( ERROR );
  }

  numByte = numPeak * sizeof( int );
  compress = (int *)malloc( numByte );
  if( compress == NULL )
  {
    fprintf( stderr, "fit_compress: error: unable to allocate memory\n" );
    free( seq );
    return( ERROR );
  }

  for( i = 0; i < numPeak; ++i )
  {
    compress[i] = 0;
  }

  numByte = numPeak * sizeof( FLOAT );
  relative_area = (FLOAT *)malloc( numByte );
  if( relative_area == NULL )
  {
    fprintf( stderr, "fit_compress: error: unable to allocate memory\n" );
    free( seq );
    free( compress );
    return( ERROR );
  }

  area = (FLOAT *)malloc( numByte );
  if( area == NULL )
  {
    fprintf( stderr, "fit_compress: error: unable to allocate memory\n" );
    free( relative_area );
    free( seq );
    free( compress );
    return( ERROR );
  }

  /*
  ** Calculate some base parameters.
  */
  baseQual = initBaseQual( first_peak, &numPeak );

  resPar7( tr_length, tr_vals, tot_vals, numPeak, baseQual );

  spcRatioBaseQual( numPeak, baseQual );

  for( i = 0; i < numPeak; ++i )
  {
    relative_area[i] = 0.0;
    area[i]          = 0.0;
  }

  /*
  ** Copy bases to sequence array.
  */
  for( i = 0, peak = first_peak; peak; peak = peak->next )
  {
    obs_peak = peak->obs_peak;
    if( obs_peak )
    {
      type = obs_peak->type;
      for( shift = 1; shift <= type; shift++ )
      {
        seq[i]           = "ACGTN"[(int)obs_peak->nuc];
        relative_area[i] = obs_peak->relative_area;
        area[i]          = obs_peak->area / (FLOAT)obs_peak->type;
        ++i;

        if( shift != type )
        {
          peak = peak->next;
        }
      }
    }
    else
    {
      seq[i]           = 'N';
      relative_area[i] = -1.0;
      area[i]          = -1.0;
      ++i;
    }
  }

  /*
  ** Calculate a base specific relative area.
  */
  numByte = numPeak * sizeof( int );
  inx = (int *)malloc( numByte );
  if( inx == NULL )
  {
    fprintf( stderr, "fit_compress: error: unable to allocate memory\n" );
    free( area );
    free( relative_area );
    free( compress );
    free( seq );
    free( baseQual );
    return( ERROR );
  }

  numByte = numPeak * sizeof( FLOAT );
  area_sum = (FLOAT *)malloc( numByte );
  if( area_sum == NULL )
  {
    fprintf( stderr, "fit_compress: error: unable to allocate memory\n" );
    free( area );
    free( relative_area );
    free( compress );
    free( seq );
    free( inx );
    free( baseQual );
    return( ERROR );
  }

  numByte = numPeak * sizeof( FLOAT );
  area_var = (FLOAT *)malloc( numByte );
  if( area_var == NULL )
  {
    fprintf( stderr, "fit_compress: error: unable to allocate memory\n" );
    free( area_sum );
    free( area );
    free( relative_area );
    free( compress );
    free( seq );
    free( inx );
    free( baseQual );
    return( ERROR );
  }

  for( i = 0; i < numPeak; ++i )
  {
    area_sum[i]      = 0.0;
    area_var[i]      = 0.0;
    relative_area[i] = 0.0;
  }

  base[0] = 'A';
  base[1] = 'C';
  base[2] = 'G';
  base[3] = 'T';
  for( j = 0; j < 4; ++j )
  {

    /*
    ** Locate peaks of base type.
    ** (Overlook Ns.)
    */
    n = 0;
    for( i = 0; i < numPeak; ++i )
    {
      if( seq[i] == base[j] )
      {
        inx[n] = i;
        if( n >= 5 )
        {
          area_sum[i] = area[i] +
                        area[inx[n-1]] +
                        area[inx[n-2]] +
                        area[inx[n-3]] +
                        area[inx[n-4]] +
                        area[inx[n-5]];
        }
        ++n;
      }
    }

    /*
    ** Calculate relative areas if there are more
    ** than six called peaks of this base type.
    */
    if( n > 12 )
    {
      for( i = 0; i < n; ++i )
      {
        if( i < 6 )
        {
          relative_area[inx[i]] = area_sum[inx[i+6]] > 0.0 ? area[inx[i]] / ( area_sum[inx[i+6]] / 6.0 ) : 0.0;
        }
        else
        if( i >= n - 6 )
        {
          relative_area[inx[i]] = area_sum[inx[i-1]] > 0.0 ? area[inx[i]] / ( area_sum[inx[i-1]] / 6.0 ) : 0.0;
        }
        else
        if( ( area_sum[inx[i-1]] + area_sum[inx[i+6]] ) > 0.0 )
        {
          relative_area[inx[i]] = area[inx[i]] / ( ( area_sum[inx[i-1]] + area_sum[inx[i+6]] ) / 12.0 );
        }
        else
        {
          /*
          ** Set relative areas to zero to remove them from consideration.
          */
          relative_area[inx[i]] = 0.0;
        }
      }
    }
    else
    {
      /*
      ** Set relative areas to zero to remove them from consideration.
      */
      for( i = 0; i < n; ++i )
      {
        relative_area[inx[i]] = 0.0;
      }
    }
  }

  free( area );
  free( inx );
  free( area_sum );

  /*
  ** Find compression motifs in sequence.
  */
  compressMotif( numPeak, seq, compress );

  /*
  ** Locate split candidates.
  */
  for( peak = first_peak, i = 0; peak; peak = peak->next, ++i )
  {
    /*
    ** Split candidate must meet the conditions
    **   1. not the same base as the preceding called base
    **      (mononucleotide run peaks tend to increase in
    **       amplitude)
    **   2. compression motif
    **   3. peak is larger than surrounding peaks
    **   4. resolution is reasonable
    **   5. peak spacing varies at least a little
    **   6. peak is not split already
    */

/*
printf( "fit_compress: %-6d  relArea: %-6.3f    resPar: %-6.3f    spcRatio: %-6.3f    type: %d    compress: %d    diff: %d\n",
        i + 1,
        relative_area[i],
        baseQual[i].resPar7,
        baseQual[i].spcRatio,
        peak->obs_peak &&
        peak->obs_peak->type,
        compress[i],
        i > 1 ? ( seq[i] != seq[i-1] ) : 0 );
*/


    if( i > 1 &&
        seq[i] != seq[i-1] &&
        compress[i] == 1 &&
        relative_area[i] > 1.6 &&
        baseQual[i].resPar7  < 0.8 &&
        baseQual[i].spcRatio > 1.1 &&
        baseQual[i].spcRatio < 1.6 &&
        peak->obs_peak &&
        peak->obs_peak->type == 1 )
    {
      /*
      ** Adjust relative area of observed peak.
      */
      fac = (float)peak->obs_peak->type / (float)( peak->obs_peak->type + 1 );
      peak->obs_peak->relative_area = fac * peak->obs_peak->relative_area;

      /*
      ** Initialize and link in new peak.
      */
      new_peak             = alloc_predicted_peak();
      init_peak( new_peak,
                 peak->pred_location,
                 peak->pred_period,
                 peak->total_signal,
                 peak->proportion_fitted );
      new_peak->prev       = peak;
      new_peak->next       = peak->next;
      peak->next           = new_peak;
      new_peak->next->prev = new_peak;
      new_peak->obs_peak   = peak->obs_peak;

      /*
      ** Adjust structure values.
      */
      peak->obs_peak->type = peak->obs_peak->type + 1;

      /*
      ** Update obs_peak->peak[MAX_NUM_SPLITS+1]
      */
      new_peak->obs_peak->peak[new_peak->obs_peak->type] = new_peak;

      /*
      ** Identify new_peak->nuc.
      */
      new_peak->nuc = new_peak->obs_peak->nuc;

      /*
      ** Move peak to inserted peak in order to keep 'i' meaningful.
      */
      peak = new_peak;
    }
  }

  /*
  ** Free memory.
  */
  free( area_var );
  free( relative_area );
  free( compress );
  free( seq );
  free( baseQual );

  return( OK );
}

