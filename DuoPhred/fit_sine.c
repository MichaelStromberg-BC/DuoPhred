/** fit_sine.c **/

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
**  Local basecalling front-end.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "phred.h"

#ifdef ANSI_C
int normalize_traces( FILE *fp, FLOAT **tr_vals, FLOAT *tot_vals,
                      int tr_length, FLOAT *scale_traces );
#else
int normalize_traces();
#endif

/*
** fit_main:
**
** local base calling routine
**
*/
#ifdef ANSI_C
int fit_main( int tr_length, FLOAT **tr_vals,
              PhredData *phredData, int *status )
#else
int fit_main( tr_length, tr_vals, phredData, status )
int tr_length;
FLOAT **tr_vals;
PhredData *phredData;
int *status;
#endif
{
  int i, j;
  int shift, type;
  int numPeak;
  int istat;
  int *output_pos;
  int *output_qual;
  int *output_qual_uncalled;
  FLOAT *tot_vals;
  FLOAT scale_traces[4];
  char *output_nucs;
  Observed_peak *obs_peak;
  Observed_peak *first_obs_peak;
  Peak *first_peak, *peak;
  Option *option;

  option = getOption();

  tot_vals = (FLOAT *)ourMalloc( tr_length * sizeof(FLOAT) );

  /*
  ** Call bases.
  */
  first_peak = fit_sine( tr_vals, tot_vals, tr_length,
                         scale_traces, &first_obs_peak,
                         phredData->compressSplitFlag,
                         phredData->chemType, &istat );
  if( istat == ERROR )
  {
    ourFree( (char *)tot_vals );
    *status = ERROR;
    return( -1 );
  }

  /*
  ** Count peaks.
  */
  numPeak = 0;
  for( i = 0, peak = first_peak; peak; peak = peak->next )
  {
    ++numPeak;
  }

  /*
  ** Allocate memory.
  */
  output_nucs = (char *)ourMalloc( numPeak * sizeof( char ) );
  output_pos  = (int *)ourMalloc( numPeak * sizeof( int ) );
  output_qual = (int *)ourMalloc( numPeak * sizeof( int ) );
  output_qual_uncalled = (int *)ourMalloc( numPeak * sizeof( int ) );

  /*
  ** Copy base calls to arrays and prepare to set base call
  ** quality values.
  */
  for( i = 0, peak = first_peak; peak; peak = peak->next )
  {
    obs_peak = peak->obs_peak;
    if( obs_peak )
    {
      type = obs_peak->type;
      for( shift = 1; shift <= type; shift++ )
      {
        peak->nuc = "ACGTN"[(int)obs_peak->nuc];
        output_nucs[i] = peak->nuc;
        output_pos[i] = obs_peak->split[type][shift];
        ++i;

        if( shift != type )
        {
          peak = peak->next;	
        }
      }
    }
    else
    {
      peak->nuc = 'N';
      output_nucs[i] = peak->nuc;
      output_pos[i]  = peak->pred_location;
      ++i;
    }
  }

  /*
  ** Set base call quality values.
  */
  setQual( tr_length, tr_vals, tot_vals,
           i, first_peak, output_qual, output_qual_uncalled, phredData,
           option->qualityValueCeilingOption,
           option->qualityValueCeiling );

  /*
  ** Load PolyData structure array if requested.
  */
  if( option->writePolyData )
  {
    loadPolyData( first_peak, scale_traces, tr_vals, &(phredData->polyData) );
  }

  /*
  ** Perform a few sanity checks.
  */
  for( j = 0; j < i; ++j )
  {
    if( output_pos[j] < 0 ||
        output_pos[j] >= tr_length )
    {
      fprintf( stderr, "fit_main: bad base position: out of range: base number: %d  pos: %d\n",
               j + 1, output_pos[j] );
      output_pos[j] = 0;
    }

    if( output_nucs[j] != 'A' &&
        output_nucs[j] != 'C' &&
        output_nucs[j] != 'G' &&
        output_nucs[j] != 'T' &&
        output_nucs[j] != 'N' )
    {
      fprintf( stderr, "fit_main: bad base call: unknown base: base number: %d  base: %c\n",
               j + 1, output_nucs[j] );
      output_nucs[j] = 'N';
    }

    if( output_qual[j] < 0 ||
        output_qual[j] > 62 )
    {
      fprintf( stderr, "fit_main: bad quality value: base number: %d  quality: %d\n", j + 1, output_qual[j] );
      output_qual[j] = 0;
    }
    
    if( output_qual_uncalled[j] < 0 ||
        output_qual_uncalled[j] > 62 )
    {
      fprintf( stderr, "fit_main: bad quality value: base number: %d  quality (uncalled): %d\n", j + 1, output_qual_uncalled[j] );
      output_qual_uncalled[j] = 0;
    }
    
  }


  /*
  ** Copy pointers.
  */
  phredData->numBase[LCL]  = i;
  phredData->base[LCL]     = output_nucs;
  phredData->baseLoc[LCL]  = output_pos;
  phredData->baseQual[LCL] = output_qual;
  phredData->baseQualUncalled[LCL] = output_qual_uncalled;

  //printf("MPS: pointers copied.\n");

  /*
  ** Free memory.
  */
  ourFree( (char *)tot_vals );
  free_obs_peak( first_obs_peak );
  free_predicted_peak();

  *status = OK;

  return( i );

}



/*
** fit_sine:
*/
#ifdef ANSI_C
Peak *fit_sine( FLOAT **tr_vals, FLOAT *tot_vals, int tr_length,
                FLOAT *scale_traces, Observed_peak **pfirst_obs_peak,
                int compressSplitFlag, int chemType, int *status )
#else
Peak *fit_sine( tr_vals, tot_vals, tr_length, scale_traces,
                pfirst_obs_peak, compressSplitFlag, chemType, status )
FLOAT **tr_vals;
FLOAT *tot_vals;
int tr_length;
FLOAT *scale_traces;
Observed_peak **pfirst_obs_peak;
int compressSplitFlag;
int chemType;
int *status;
#endif
{
  int npk;
  int i, j, n_fixed;
  int begGoodTrace;
  int endGoodTrace;
  FLOAT **sig_func;
  Peak *first_peak;
  Observed_peak *first_obs_peak;
  LocPeak *locPeak;
  FILE *fp;
  Option *option;

  option = getOption();

  fp = stdout; 

  sig_func = (FLOAT **)ourMalloc( 4 * sizeof( FLOAT * ) );
  for( i = 0; i < 4; ++i )
  {
    sig_func[i] = (FLOAT *)ourMalloc( tr_length * sizeof( FLOAT ) );
    for( j = 0; j < tr_length; ++j )
    {
      sig_func[i][j] = tr_vals[i][j];
    }
  }

  /*
  ** Normalize the four traces relative to each other.
  */
  if( option->normalize && tr_length > 2500 )
  {
    normalize_traces( fp, tr_vals, tot_vals, tr_length, scale_traces );
  }
  else
  {
    /*
    ** The polyData modules need informative scaling values.
    */
    for( j = 0; j < 4; ++j )
    {
      scale_traces[j] = -1.0;
    }
  }

  /*
  ** Define a "sky-line" projection of the four traces.
  */
  makeTotVals( tr_vals, tot_vals, tr_length );

  /*
  ** Find the maximum trace value of each of the four traces.
  */
  setMaxVal( tr_length, tr_vals );

  /*
  ** Evaluate the trace quality before calling bases.
  */
  locPeak = evaluateTrace( tr_length, tr_vals, tot_vals, 200,
                           &npk, &begGoodTrace, &endGoodTrace );

  if( locPeak == NULL )
  {
    for( i = 0; i < 4; ++i )
    {
      ourFree( (char *)sig_func[i] );
    }
    ourFree( (char *)sig_func );
    *status = ERROR;
    return( NULL );
  }

  /*
  ** Locate where peaks are expected.
  */
  first_peak = make_predicted_peaks( fp, tr_vals, tr_length, tot_vals,
                                     npk, locPeak,
                                     option->beginPeakPredictionOption,
                                     option->beginPeakPredictionPoint );


  /*
  ** Locate observed peaks (trace peaks).
  */
  first_obs_peak = find_observed_peaks( tr_vals, tr_length );
  *pfirst_obs_peak = first_obs_peak;

  /*
  ** Perform base-calling.
  */
  n_fixed = fit_peaks( tr_length, tr_vals, first_peak, first_obs_peak );

  /*
  ** Try to fill in obvious gaps in the base calling.
  */
  capture_missed_peaks1( tr_length, first_peak, first_obs_peak, sig_func, begGoodTrace, endGoodTrace );
  capture_missed_peaks2( tr_length, tr_vals, tot_vals, first_peak, first_obs_peak );

  if( compressSplitFlag == 1 && chemType == PRIMER_CHEM )
  {
    /*
    ** Find compressions.
    */
    fit_compress( first_peak, tr_length, tr_vals, tot_vals );
  }

  /*
  ** Free memory.
  */
  for( i = 0; i < 4; ++i )
  {
    ourFree( (char *)sig_func[i] );
  }
  ourFree( (char *)sig_func );
  ourFree( (char *)locPeak );

  *status = OK;

  return( first_peak );
}  


