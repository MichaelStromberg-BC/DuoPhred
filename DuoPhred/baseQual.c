/** baseQual.c **/

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
** Base quality functions.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "phred.h"


/*
** Allocate and initialize elements of
** base quality structure array.
*/
#ifdef ANSI_C
BaseQual *initBaseQual( Peak *first_peak, int *numBase )
#else
BaseQual *initBaseQual( first_peak, numBase )
Peak *first_peak;
int *numBase;
#endif
{
  int i;
  int type;
  int shift;
  Peak *peak;
  Observed_peak *obs_peak;
  BaseQual *baseQual;

  i = 0;
  for( peak = first_peak; peak; peak = peak->next )
  {
    ++i;
  }

  if( ( baseQual = (BaseQual *)ourMalloc( i * sizeof( BaseQual ) ) ) == NULL )
  {
    return( NULL );
  }

  i = 0;
  for( peak = first_peak; peak; peak = peak->next )
  {
    obs_peak = peak->obs_peak;
    if( obs_peak )
    {
      type = obs_peak->type;
      for( shift = 1; shift <= type; ++shift )
      {
        baseQual[i].nuc          = obs_peak->nuc;
        baseQual[i].loc          = obs_peak->split[type][shift];
        baseQual[i].type         = type;
        baseQual[i].relativeArea = obs_peak->relative_area;
        baseQual[i].peak         = peak;
        if( peak->best_uncalled_peak )
        {
          baseQual[i].locUncalled = peak->best_uncalled_peak->split[1][1];
          baseQual[i].nucUncalled = peak->best_uncalled_peak->nuc;
        }
        else
        {
          baseQual[i].locUncalled = -1;
          baseQual[i].nucUncalled = -1;
        }
        ++i;
        if( shift != type )
        {
          peak = peak->next;
        }
      }
    }
    else
    {
      baseQual[i].nuc          = 4;
      baseQual[i].loc          = peak->pred_location;
      baseQual[i].type         = 0;
      baseQual[i].relativeArea = 0.0;
      baseQual[i].peak         = peak;
      if( peak->best_uncalled_peak )
      {
        baseQual[i].locUncalled  = peak->best_uncalled_peak->split[1][1];
        baseQual[i].nucUncalled  = peak->best_uncalled_peak->nuc;
      }
      else
      {
        baseQual[i].locUncalled  = -1;
        baseQual[i].nucUncalled  = -1;
      }
      ++i;
    }
  }

  *numBase = i;

  return( baseQual );

}




/*
** Set maxDown elements in base quality structure array.
*/
#ifdef ANSI_C
int maxDownBaseQual( int tr_length, FLOAT **tr_vals, FLOAT *tot_vals, int numBase, BaseQual *baseQual )
#else
int maxDownBaseQual( tr_length, tr_vals, tot_vals, numBase, baseQual )
int tr_length;
FLOAT **tr_vals;
FLOAT *tot_vals;
int numBase;
BaseQual *baseQual;
#endif
{
  int i, j, k;
  int itmp;
  int loc1, nuc1;
  int loc2, nuc2;
  FLOAT minPair;

  for( i = 0; i < ( numBase - 1 ); ++i )
  {
    /*
    ** Pair of adjacent bases.
    */
    loc1 = baseQual[i].loc < tr_length ? baseQual[i].loc : tr_length - 1;
    nuc1 = baseQual[i].nuc;
    loc2 = baseQual[i+1].loc < tr_length ? baseQual[i+1].loc : tr_length - 1;
    nuc2 = baseQual[i+1].nuc;

    if( loc1 > loc2 )
    {
      itmp = loc1;
      loc1 = loc2;
      loc2 = itmp;

      itmp = nuc1;
      nuc1 = nuc2;
      nuc2 = itmp;
    }

    /*
    ** Find smallest of these two base peaks.
    */
    if( nuc1 < 4 && nuc2 < 4 )
    {
      minPair = ( tr_vals[nuc1][loc1] < tr_vals[nuc2][loc2] ) ?
                tr_vals[nuc1][loc1] : tr_vals[nuc2][loc2];
    }
    else
    if( nuc1 == 4 && nuc2 == 4 )
    {
      minPair = -1.0;
    }
    else
    if( nuc1 == 4 )
    {
      minPair = tr_vals[nuc2][loc2];
    }
    else
    if( nuc2 == 4 )
    {
      minPair = tr_vals[nuc1][loc1];
    }
    else
    {
      printf( "unrecognized base type\n" );
      minPair = 0.0;
    }

    /*
    ** Find lowest trace value between two peaks.
    */
    k = loc1;
    for( j = loc1; j <= loc2; ++j )
    {
      if( tot_vals[j] < tot_vals[k] )
      {
        k = j;
      }
    }
    baseQual[i].maxDownF1 = tot_vals[k];

    /*
    ** Ratio of lowest value-to-smallest peak.
    */
    if( minPair > 0.0 )
    {
      baseQual[i].maxDownF1 /= minPair;
    }
    else
    {
      baseQual[i].maxDownF1 = -1.0;
    }
  }
  baseQual[numBase-1].maxDownF1 = baseQual[numBase-2].maxDownF1;

  /*
  ** Count number of bases on either side of base of interest
  ** that have dips in the trace between peaks.
  */
  for( i = 0; i < numBase; ++i )
  {
    j = 1;
    if( baseQual[i].relativeArea > 0.0 &&
        baseQual[i].type == 1 &&
        baseQual[i].nuc != 'N' )
    {
      while( 1 )
      {
        if( i - j < 0 ||
            baseQual[i-j].maxDownF1 > 0.99 ||
            baseQual[i-j].relativeArea <= 0.0 ||
            baseQual[i-j].type != 1 ||
            baseQual[i-j].nuc == 'N' )
        {
          break;
        }
        if( i + j >= numBase ||
            baseQual[i+j-1].maxDownF1 > 0.99 ||
            baseQual[i+j].relativeArea <= 0.0 ||
            baseQual[i+j].type != 1 ||
            baseQual[i+j].nuc == 'N' )
        {
          break;
        }
        ++j;
      }
    }
    baseQual[i].maxDown = (float)( -( j - 1 ) );
  }

  return( OK );
}


/*
** Set spcRatio elements in base quality structure array.
*/
#ifdef ANSI_C
int spcRatioBaseQual( int numBase, BaseQual *baseQual )
#else
int spcRatioBaseQual( numBase, baseQual )
int numBase;
BaseQual *baseQual;
#endif
{
  int i, j;
  FLOAT min, max;

  /*
  ** Calculate space between pairs of peak locations.
  */
  for( i = 0; i < ( numBase - 1 ); ++i )
  {
    baseQual[i].space = baseQual[i+1].loc - baseQual[i].loc;
  }
  baseQual[numBase-1].space = baseQual[numBase-2].space;

  /*
  ** Calculate space ratio.
  */
  for( i = 3; i < ( numBase - 2 ); ++i )
  {
    min = baseQual[i-3].space;
    max = baseQual[i-3].space;
    baseQual[i].basMaxSpc = i - 3;
    for( j = -2; j <= 2; ++j )
    {
      if( baseQual[i+j].space < min )
      {
        min = baseQual[i+j].space;
      }
      if( baseQual[i+j].space > max )
      {
        max = baseQual[i+j].space;
        baseQual[i+j].basMaxSpc = i + j;
      }
    }
    baseQual[i].spcRatio = min > 0.0 ? max / min : 100.0;
  }

  /*
  ** Extend values to first and last three array elements.
  */
  for( i = 0; i < 3; ++i )
  {
    baseQual[i].spcRatio = baseQual[3].spcRatio;
  }
  for( i = numBase - 2; i <= ( numBase - 1 ); ++i )
  {
    baseQual[i].spcRatio = baseQual[numBase-3].spcRatio;
  }

  return( OK );
}


/*
** Set maxRatio3 and maxRatio7 elements in base quality
** structure array.
*/
#ifdef ANSI_C
int maxRatioBaseQual( int tr_length, FLOAT **tr_vals, FLOAT **otr_vals, int un_norm, int numBase, BaseQual *baseQual )
#else
int maxRatioBaseQual( tr_length, tr_vals, otr_vals, un_norm, numBase, baseQual )
int tr_length;
FLOAT **tr_vals;
FLOAT **otr_vals;
int un_norm;
int numBase;
BaseQual *baseQual;
#endif
{
  int i, j, k;
  int loc, nuc;
  int numBytes;
  FLOAT min, omin;
  FLOAT max, omax;
  FLOAT maxRatio, omaxRatio;
  FLOAT mpsMaxRatio; 
	//FLOAT* maxRatios3;
	//FLOAT* maxRatios7;
	FLOAT difference;
	/*
	//
  // Allocate and initialize memory.
  //
  numBytes = numBase * sizeof(FLOAT);
  
  maxRatios3 = (FLOAT *)malloc(numBytes);
  if(maxRatios3 == NULL) {
    fprintf( stderr, "maxRatios3: error: unable to allocate memory\n" );
    return( ERROR );
  }
	
	maxRatios7 = (FLOAT *)malloc(numBytes);
  if(maxRatios7 == NULL) {
    fprintf( stderr, "maxRatios7: error: unable to allocate memory\n" );
    return( ERROR );
  }
	
	printf("MPS: maxRatioBaseQual called.\n");
	*/
  /*
  ** Find the amplitude of the called peaks.
  */
  for( i = 0; i < numBase; ++i )
  {
    if( baseQual[i].nuc != 4 )
    {
      baseQual[i].maxCalled  = tr_vals[baseQual[i].nuc][baseQual[i].loc];
      baseQual[i].omaxCalled = otr_vals[baseQual[i].nuc][baseQual[i].loc];
    }
    else
    {
      baseQual[i].maxCalled  = 0.0;
      baseQual[i].omaxCalled = 0.0;
    }
  }

  /*
  ** Find the amplitude of the uncalled peak, or, if there is no
  ** uncalled peak, the highest amplitude trace that is not the
  ** trace of the called base at the location of the called base.
  */
  for( i = 0; i < numBase; ++i )
  {
    if( baseQual[i].locUncalled >= 0 )
    {
      nuc = baseQual[i].nucUncalled;
      loc = baseQual[i].locUncalled;
      if( nuc >= 0 && nuc < 4 )
      {
        baseQual[i].maxUncalled  = tr_vals[nuc][loc];
        baseQual[i].omaxUncalled = otr_vals[nuc][loc];
        

      }
      else
      {
        baseQual[i].maxUncalled  = -1.0;
        baseQual[i].omaxUncalled = -1.0;
        printf( "unexpected peak type: %d\n", nuc );
      }
    }
    else
    {
      if( baseQual[i].nuc == 4 )
      {
        baseQual[i].maxUncalled  = -1.0;
        baseQual[i].omaxUncalled = -1.0;
        continue;
      }

      loc = baseQual[i].loc;
      k = ( baseQual[i].nuc != 0 ) ? 0 : 1;
      for( j = 0; j < 4; ++j )
      {
        if( j != baseQual[i].nuc &&
            tr_vals[j][loc] > tr_vals[k][loc] )
        {
          k = j;
        }
      }
      baseQual[i].maxUncalled = tr_vals[k][loc];

      loc = baseQual[i].loc;
      k = ( baseQual[i].nuc != 0 ) ? 0 : 1;
      for( j = 0; j < 4; ++j )
      {
        if( j != baseQual[i].nuc &&
            otr_vals[j][loc] > otr_vals[k][loc] )
        {
          k = j;
        }
      }
      baseQual[i].omaxUncalled = otr_vals[k][loc];
    }
  }
    
  /*
  ** Calculate the amplitude ratio for three peaks.
  */
  for( i = 0; i < ( numBase - 2 ); ++i )
  {
    min  = baseQual[i].maxCalled;
    max  = baseQual[i].maxUncalled;
    omin = baseQual[i].omaxCalled;
    omax = baseQual[i].omaxUncalled;

		//printf("MPS: base %d min %f max %f omin %f omax %f\n",i,min,max,omin,omax);
		
    for( j = 1; j < 3; ++j )
    {
      if( baseQual[i+j].maxCalled < min )
      {
        min = baseQual[i+j].maxCalled;
      }
      if( baseQual[i+j].maxUncalled > max )
      {
        max = baseQual[i+j].maxUncalled;
      }

      if( baseQual[i+j].omaxCalled < omin )
      {
        omin = baseQual[i+j].omaxCalled;
      }
      if( baseQual[i+j].omaxUncalled > omax )
      {
        omax = baseQual[i+j].omaxUncalled;
      }

    }
    
    maxRatio  = min  > 0.0 ? max  / min  : 100.0;
    
    if( un_norm )
    {
      omaxRatio = omin > 0.0 ? omax / omin : 100.0;
      baseQual[i+1].maxRatio3 = maxRatio > omaxRatio ? maxRatio : omaxRatio;
    }
    else
    {
    	// WE USUALLY DO THIS
      baseQual[i+1].maxRatio3 = maxRatio;
    }
  }
  baseQual[0].maxRatio3 = baseQual[1].maxRatio3;
  baseQual[numBase-1].maxRatio3 = baseQual[numBase-2].maxRatio3;

  /*
  ** Calculate the amplitude ratio for seven peaks.
  */
  for( i = 0; i < ( numBase - 6 ); ++i )
  {
    min  = baseQual[i].maxCalled;
    max  = baseQual[i].maxUncalled;
    omin = baseQual[i].omaxCalled;
    omax = baseQual[i].omaxUncalled;
    for( j = 1; j < 7; ++j )
    {
      if( baseQual[i+j].maxCalled < min )
      {
        min = baseQual[i+j].maxCalled;
      }
      if( baseQual[i+j].maxUncalled > max )
      {
        max = baseQual[i+j].maxUncalled;
      }

      if( baseQual[i+j].omaxCalled < omin )
      {
        omin = baseQual[i+j].omaxCalled;
      }
      if( baseQual[i+j].omaxUncalled > omax )
      {
        omax = baseQual[i+j].omaxUncalled;
      }
    }
    maxRatio  = min  > 0.0 ? max  / min  : 100.0;
    if( un_norm )
    {
      omaxRatio = omin > 0.0 ? omax / omin : 100.0;
      baseQual[i+3].maxRatio7 = maxRatio > omaxRatio ? maxRatio : omaxRatio;
    }
    else
    {
      baseQual[i+3].maxRatio7 = maxRatio;
    }
  }
  baseQual[0].maxRatio7 = baseQual[3].maxRatio7;
  baseQual[1].maxRatio7 = baseQual[3].maxRatio7;
  baseQual[2].maxRatio7 = baseQual[3].maxRatio7;
  baseQual[numBase-3].maxRatio7 = baseQual[numBase-4].maxRatio7;
  baseQual[numBase-2].maxRatio7 = baseQual[numBase-4].maxRatio7;
  baseQual[numBase-1].maxRatio7 = baseQual[numBase-4].maxRatio7;

	//
	// build my maximum ratios for 3 peaks using the uncalled
	// nucleotide
	//
	for(i = 0; i < (numBase - 2); i++) {
		baseQual[i+1].maxRatio3Uncalled = getMaxRatioNU(i, 3, tr_vals, baseQual);
	}
	
	baseQual[0].maxRatio3Uncalled = baseQual[1].maxRatio3Uncalled;
	baseQual[numBase-1].maxRatio3Uncalled = baseQual[numBase-2].maxRatio3Uncalled;

	//
	// build my maximum ratios for 7 peaks using the uncalled nucleotide
	//
  
  for(i = 0; i < (numBase - 6); i++) {
  	baseQual[i+3].maxRatio7Uncalled = getMaxRatioNU(i, 7, tr_vals, baseQual);
		//maxRatios7[i+3] = getMaxRatioN(i, 7, tr_vals, baseQual);
	}
	
	baseQual[0].maxRatio7Uncalled = baseQual[3].maxRatio7Uncalled;
	baseQual[1].maxRatio7Uncalled = baseQual[3].maxRatio7Uncalled;
	baseQual[2].maxRatio7Uncalled = baseQual[3].maxRatio7Uncalled;
	baseQual[numBase-3].maxRatio7Uncalled = baseQual[numBase-4].maxRatio7Uncalled;
	baseQual[numBase-3].maxRatio7Uncalled = baseQual[numBase-4].maxRatio7Uncalled;
	baseQual[numBase-1].maxRatio7Uncalled = baseQual[numBase-4].maxRatio7Uncalled;

  return( OK );
}

// MPS: compare function
int compare(void *context, const void *v1, const void *v2) {
		return (*(FLOAT*)v1 - *(FLOAT *)v2);
}

// MPS: function to calculate the maximum ratio for the next
//      N peaks.
FLOAT getMaxRatioN(int loc, int windowLen, FLOAT **tr_vals, BaseQual *baseQual) {
	
	int i, j;
	int numBytes;
	int locCalled, locUncalled;
	int nucCalled, nucUncalled;
	int maxNuc;
	FLOAT min, max;
	
	FLOAT* calledTraceVals;
	FLOAT* uncalledTraceVals;
	
	FLOAT maxRatio;
	
	//
	// allocate and initialize memory	
  //
  numBytes = windowLen * sizeof(FLOAT);
  
  calledTraceVals = (FLOAT *)malloc(numBytes);
  if(calledTraceVals == NULL) {
    fprintf(stderr, "ERROR in getMaxRatioN: unable to allocate memory for calledTraceVals.\n" );
    return(ERROR);
  }
  
  uncalledTraceVals = (FLOAT *)malloc(numBytes);
  if(uncalledTraceVals == NULL) {
    fprintf(stderr, "ERROR in getMaxRatioN: unable to allocate memory for uncalledTraceVals.\n" );
    return(ERROR);
  }
	
	// get the minimum called values and maximum uncalled values 
	// for the next N bases
	for(i = 0; i < windowLen; i++) {
		
		// set our location
		locCalled   = baseQual[loc+i].loc;
		locUncalled = baseQual[loc+i].locUncalled;
		
		// set our bases
		nucCalled   = baseQual[loc+i].nuc;
		nucUncalled = baseQual[loc+i].nucUncalled;
		
		// if we don't have an uncalled peak, use the main peak location
		if(locUncalled < 0) locUncalled = locCalled;
		
		if(nucUncalled < 0) {
			
			// get the highest non-called peak
			maxNuc = -1;
			max    = -1;
			for(j = 0; j < 4; j++) {
				
				// skip the called nucleotide
				if(j == nucCalled) continue;
				
				if(tr_vals[j][locUncalled] > max) {
					max = tr_vals[j][locUncalled];
					maxNuc = j;
				}
			}

			// our new uncalled nucleotide			
			nucUncalled = maxNuc;
		}
				
		// set our values
		calledTraceVals[i]   = tr_vals[nucCalled][locCalled];
		uncalledTraceVals[i] = tr_vals[nucUncalled][locUncalled];
	}
	
	// do some diagnostics
	//printf("before sort: ctv0: %f, ctv1: %f, ctv2: %f\n",calledTraceVals[0],calledTraceVals[1],calledTraceVals[2]);
	//printf("             utv0: %f, utv1: %f, utv2: %f\n",uncalledTraceVals[0],uncalledTraceVals[1],uncalledTraceVals[2]);
	
	// do some quick sorts
	qsort_s(calledTraceVals, windowLen, sizeof(FLOAT), compare, NULL);	
	qsort_s(uncalledTraceVals, windowLen, sizeof(FLOAT), compare, NULL);
	
	// do some diagnostics
	//printf("after sort: ctv0: %f, ctv1: %f, ctv2: %f\n",calledTraceVals[0],calledTraceVals[1],calledTraceVals[2]);
	//printf("            utv0: %f, utv1: %f, utv2: %f\n",uncalledTraceVals[0],uncalledTraceVals[1],uncalledTraceVals[2]);
	
	// get the minimum of the called values, and the maximum of the uncalled values
	min  = calledTraceVals[0];
	max  = uncalledTraceVals[windowLen-1];
	
	//printf("# compared: %d, min: %f, max: %f\n\n",numCompared,min,max);
	
	// calculate the maximum ratio for N peaks
	maxRatio = min > 0.0 ? max  / min : 100.0;
	
	// free memory
	if(calledTraceVals) free(calledTraceVals);
	if(uncalledTraceVals) free(uncalledTraceVals);
	
	// return the maximum ratio
	return maxRatio;
}

// MPS: function to calculate the maximum ratio for the next
//      N peaks. Modified to handle an uncalled primary base
//
//      The trick here is that the middle sample needs to reference
//      the uncalled base instead of the called base.
FLOAT getMaxRatioNU(int loc, int windowLen, FLOAT **tr_vals, BaseQual *baseQual) {
	
	int i, j;
	int numBytes;
	int locCalled, locUncalled;
	int nucCalled, nucUncalled;
	int maxNuc;
	int pivotPoint;
	FLOAT min, max;
	
	FLOAT* calledTraceVals;
	FLOAT* uncalledTraceVals;
	
	FLOAT maxRatio;
	
	//
	// allocate and initialize memory	
  //
  numBytes = windowLen * sizeof(FLOAT);
  
  calledTraceVals = (FLOAT *)malloc(numBytes);
  if(calledTraceVals == NULL) {
    fprintf(stderr, "ERROR in getMaxRatioN: unable to allocate memory for calledTraceVals.\n" );
    return(ERROR);
  }
  
  uncalledTraceVals = (FLOAT *)malloc(numBytes);
  if(uncalledTraceVals == NULL) {
    fprintf(stderr, "ERROR in getMaxRatioN: unable to allocate memory for uncalledTraceVals.\n" );
    return(ERROR);
  }
  
  // set our pivot point (should be the middle sample in an odd window length)
  pivotPoint = (int)(windowLen / 2);
  //printf("MPS: pivot point: %d\n",pivotPoint);
	
	// get the minimum called values and maximum uncalled values 
	// for the next N bases
	for(i = 0; i < windowLen; i++) {
		
		// set our location
		locCalled   = baseQual[loc+i].loc;
		locUncalled = baseQual[loc+i].locUncalled;
		
		// set our bases
		nucCalled   = baseQual[loc+i].nuc;
		nucUncalled = baseQual[loc+i].nucUncalled;
		
		// if we don't have an uncalled peak, use the main peak location
		if(locUncalled < 0) locUncalled = locCalled;
		
		if(nucUncalled < 0) {
			
			// get the highest non-called peak
			maxNuc = -1;
			max    = -1;
			for(j = 0; j < 4; j++) {
				
				// skip the called nucleotide
				if(j == nucCalled) continue;
				
				if(tr_vals[j][locUncalled] > max) {
					max = tr_vals[j][locUncalled];
					maxNuc = j;
				}
			}

			// our new uncalled nucleotide			
			nucUncalled = maxNuc;
		}
		
		// NOTE: so at this point we have defined locCalled, locUncalled, nucCalled, and nucUncalled
		
		// if this is not the uncalled pivot point
		if(i != pivotPoint) {
		
			// set our values
			// hmmm sometimes we have nucleotide 4

			if(nucCalled == 4) {
			
				calledTraceVals[i]   = 0.5;
				uncalledTraceVals[i] = 0.6;
			
			} else {

				calledTraceVals[i]   = tr_vals[nucCalled][locCalled];
				uncalledTraceVals[i] = tr_vals[nucUncalled][locUncalled];
			}
			//printf("- ctv: %f, utv: %f\n\n",calledTraceVals[i],uncalledTraceVals[i]);
		
		// this is the uncalled pivot point
		} else {
			
			// we need to find the second highest uncalled peak
			maxNuc = -1;
			max    = -1;
	
			for(j = 0; j < 4; j++) {
		
				// don't add the called nucleotide
				if(j == nucCalled) continue;
			
				// don't add the uncalled nucleotide either
				if(j == nucUncalled) continue;
			
				if(tr_vals[j][locUncalled] > max) {
					max = tr_vals[j][locUncalled];
					maxNuc = j;
				}
			}
			
			/*
			printf("called base: %c, uncalled base: %c\n","ACGTN"[nucCalled],"ACGTN"[nucUncalled]);
			printf("tv0: %f, tv1: %f, tv2: %f, tv3: %f\n",tr_vals[0][locUncalled],tr_vals[1][locUncalled],tr_vals[2][locUncalled],tr_vals[3][locUncalled]);
			printf("secondary uncalled base: %c - %f\n","ACGTN"[maxNuc],max);
			*/
			
			// set our values
			calledTraceVals[i]   = tr_vals[nucUncalled][locUncalled];
			uncalledTraceVals[i] = tr_vals[maxNuc][locUncalled];
			
			//printf("- ctv: %f, utv: %f\n\n",calledTraceVals[i],uncalledTraceVals[i]);
			
		}
	}
	
	/*
	// do some diagnostics
	for(i = 0; i < windowLen; i++) {
		printf("ctv%d: %f ",i,calledTraceVals[i]);
	}
	printf("\n");
	
	for(i = 0; i < windowLen; i++) {
		printf("utv%d: %f ",i,uncalledTraceVals[i]);
	}
	printf("\n");
	*/
	
	// do some quick sorts
	qsort_s(calledTraceVals, windowLen, sizeof(FLOAT), compare, NULL);	
	qsort_s(uncalledTraceVals, windowLen, sizeof(FLOAT), compare, NULL);
		
	// get the minimum of the called values, and the maximum of the uncalled values
	min  = calledTraceVals[0];
	max  = uncalledTraceVals[windowLen-1];
	
	//printf("# compared: %d, min: %f, max: %f\n\n",numCompared,min,max);
	
	// calculate the maximum ratio for N peaks
	maxRatio = min > 0.0 ? max  / min : 100.0;
	
	// free memory
	if(calledTraceVals) free(calledTraceVals);
	if(uncalledTraceVals) free(uncalledTraceVals);
	
	// return the maximum ratio
	return maxRatio;
}















// MPS: function to evaluate a max 3 peak for an uncalled peak
FLOAT getMaxRatio3U(int loc, FLOAT **tr_vals, FLOAT **otr_vals, BaseQual *baseQual) {
	
	int i, nucCounter;
	int locCalled, locUncalled;
	int nucCalled, nucUncalled;
	
	FLOAT min, max, omin, omax;
	
	FLOAT calledTraceVals[3];
	FLOAT uncalledTraceVals[3];
	FLOAT oCalledTraceVals[3];
	FLOAT oUncalledTraceVals[3];
	
	FLOAT maxTraceVal;
	int maxNuc;
	FLOAT maxRatio3U;
	
	// set our location
	locCalled   = baseQual[loc-1].loc;
	locUncalled = baseQual[loc-1].locUncalled;
		
	// set our bases
	nucCalled   = baseQual[loc-1].nuc;
	nucUncalled = baseQual[loc-1].nucUncalled;
	
	// So in this case we want to add the tr_vals[nucUncalled][locUncalled] to the calledTraceVals
	// we also want to find the largest uncalled trace at the locUncalled
	
	// find the nucleotide with the highest trace value at locUncalled that is neither the called
	// base (nucCalled) nor the uncalled base (nucUncalled)
	nucCounter = 0;
	maxTraceVal = 0;
	maxNuc = -1;
	
	for(i=0;i<3;i++) {
		
		// don't add the called nucleotide
		if(i == nucCalled) continue;
			
		// don't add the uncalled nucleotide either
		if(i == nucUncalled) continue;
			
		if(tr_vals[i][locUncalled] > maxTraceVal) {
			maxTraceVal = tr_vals[i][locUncalled];
			maxNuc = i;
		}
	}
	
	// at this point, treat nucUncalled like nucCalled and maxNuc like nucUncalled
	calledTraceVals[0]   = tr_vals[nucUncalled][locUncalled];
	uncalledTraceVals[0] = tr_vals[maxNuc][locUncalled];
			
	for(i = 1; i < 3; i++) {
		
		// set our location
		locCalled   = baseQual[loc+i].loc;
		locUncalled = baseQual[loc+i].locUncalled;
		
		// set our bases
		nucCalled   = baseQual[loc+i].nuc;
		nucUncalled = baseQual[loc+i].nucUncalled;
		
		// set our values
		calledTraceVals[i]  = tr_vals[nucCalled][locCalled];
		uncalledTraceVals[i]  = tr_vals[nucUncalled][locUncalled];
	}
	
	// do some quick sorts
	qsort_s(calledTraceVals, 3, sizeof(FLOAT), compare, NULL);	
	qsort_s(uncalledTraceVals, 3, sizeof(FLOAT), compare, NULL);
	
	//printf("ctv0: %f, ctv1: %f, ctv2: %f\n",calledTraceVals[0],calledTraceVals[1],calledTraceVals[2]);
	//printf("utv0: %f, utv1: %f, utv2: %f\n",uncalledTraceVals[0],uncalledTraceVals[1],uncalledTraceVals[2]);

	// get the minimum of the called values, and the maximum of the uncalled values
	min  = calledTraceVals[0];
	max  = uncalledTraceVals[2];
	
	// calculate the maxRatio3
	maxRatio3U = min > 0.0 ? max  / min : 100.0;
		
	return maxRatio3U;
}

#define FMAX(a,b)	( ( (a) > (b) ) ? (a) : (b) )

#ifdef ANSI_C
int compressParameter( int numBase, BaseQual *baseQual )
#else
int compressParameter( numBase, baseQual )
int numBase;
BaseQual *baseQual;
#endif
{
  int i;
  int numByte;
  FLOAT max1, max2;
  int *compress;
  char *seq;

  /*
  ** Allocate and initialize memory.
  */
  numByte = numBase * sizeof( int );
  compress = (int *)malloc( numByte );
  if( compress == NULL )
  {
    fprintf( stderr, "compressParameter: error: unable to allocate memory\n" );
    return( ERROR );
  }
  numByte = numBase * sizeof( char );
  seq = (char *)malloc( numByte );
  if( seq == NULL )
  {
    fprintf( stderr, "compressParameter: error: unable to allocate memory\n" );
    free( compress );
    return( ERROR );
  }
  for( i = 0; i < numBase; ++i )
  {
    seq[i] = "ACGTN"[baseQual[i].nuc];
    compress[i] = 0;
  }


  if( compressMotif( numBase, seq, compress ) == ERROR )
  {
    free( seq );
    free( compress );
    return( ERROR );
  }

  /*
  ** Calculate max relative area in window of three peaks.
  */
  max1 = FMAX( baseQual[0].relativeArea, baseQual[1].relativeArea );
  max1 = max1 < 500.0 ? max1 : 500.0;
  baseQual[0].compMaxArea = max1;
  for( i = 1; i < numBase - 1; ++i )
  {
    max2 = FMAX( baseQual[i].relativeArea, baseQual[i+1].relativeArea );
    max2 = max2 < 500.0 ? max2 : 500.0;
    baseQual[i].compMaxArea = FMAX( max1, max2 );
    max1 = max2;
  }
  baseQual[numBase-1].compMaxArea = max2;

  /*
  ** Set compPar.
  */
  for( i = 0; i < numBase; ++i )
  {
    if( baseQual[i].nuc != 'N' )
    {
      if( compress[i] == 1 ||
          ( i > 0 && compress[i-1] == 1 ) ||
          ( i < ( numBase - 1 ) && compress[i+1] == 1 ) )
      {
        baseQual[i].compPar = baseQual[i].compMaxArea;
      }
      else
      {
        baseQual[i].compPar = 0.0;
      }
    }
    else
    {
      baseQual[i].compPar = 500.0;
    }
  }

  /*
  ** Free memory.
  */
  free( seq );
  free( compress );

  return( OK );
}



#ifdef ANSI_C
int resPar7( int tr_length, FLOAT **tr_vals, FLOAT *tot_vals,
             int numBase, BaseQual *baseQual )
#else
int resPar7( tr_length, tr_vals, tot_vals, numBase, baseQual )
int tr_length;
FLOAT **tr_vals;
FLOAT *tot_vals;
int numBase;
BaseQual *baseQual;
#endif
{
  int i, j;
  int numByte;
  int lo, hi;
  int win_wid;
  int win_edg;
  int notN;
  FLOAT min;
  FLOAT min_peak;
  FLOAT max_valley;
  FLOAT *valley;

  /*
  ** Find cases in which the uncalled peak is the largest peak.
  */
  for( i = 0; i < numBase; ++i )
  {
    if( baseQual[i].locUncalled >= 0 &&
        baseQual[i].nucUncalled >= 0 &&
        baseQual[i].nucUncalled < 4 &&
        ( tot_vals[baseQual[i].locUncalled] -
          tr_vals[baseQual[i].nucUncalled][baseQual[i].locUncalled] ) < 0.01 )
    {
      baseQual[i].flgUncalled = 1;
    }
    else
    {
      baseQual[i].flgUncalled = 0;
    }
  }

  /*
  ** Allocate memory for minima between peaks.
  */
  numByte = numBase * sizeof( FLOAT );
  valley  = (FLOAT *)malloc( numByte );
  if( valley == NULL )
  {
    fprintf( stderr, "resPar7: error: unable to allocate memory\n" );
    return( ERROR );
  }

  /*
  ** Find minimum between peaks.
  */
  for( i = 0; i < ( numBase - 1 ); ++i )
  {
    min = tot_vals[baseQual[i].loc];
    for( j = baseQual[i].loc + 1;
         j <= baseQual[i+1].loc;
         ++j )
    {
      if( tot_vals[j] < min )
      {
        min = tot_vals[j];
      }
    }
    valley[i] = min;
  }
  valley[numBase-1] = valley[numBase-2];

  /*
  ** Account for possible uncalled peak between peaks.
  */
  for( i = 0; i < ( numBase - 1 ); ++i )
  {
    /*
    ** Check for uncalled peak.
    */
    if( baseQual[i].locUncalled >= 0 &&
        baseQual[i+1].locUncalled >= 0 )
    {
      if( baseQual[i].locUncalled > baseQual[i].loc &&
          baseQual[i].locUncalled < baseQual[i+1].loc &&
          baseQual[i].flgUncalled &&
          baseQual[i+1].locUncalled > baseQual[i].loc &&
          baseQual[i+1].locUncalled < baseQual[i+1].loc &&
          baseQual[i+1].flgUncalled )

      {
        if( tot_vals[baseQual[i].locUncalled] >
            tot_vals[baseQual[i+1].locUncalled] )
        {
          valley[i] = tot_vals[baseQual[i].locUncalled];
        }
        else
        {
          valley[i] = tot_vals[baseQual[i+1].locUncalled];
        }
      }
      else
      if( baseQual[i].locUncalled > baseQual[i].loc &&
          baseQual[i].locUncalled < baseQual[i+1].loc &&
          baseQual[i].flgUncalled )
      {
        valley[i] = tot_vals[baseQual[i].locUncalled];
      }
      else
      if( baseQual[i+1].locUncalled > baseQual[i].loc &&
          baseQual[i+1].locUncalled < baseQual[i+1].loc &&
          baseQual[i+1].flgUncalled )
      {
        valley[i] = tot_vals[baseQual[i+1].locUncalled];
      }
    }
    else
    if( baseQual[i].locUncalled >= 0 &&
        baseQual[i].locUncalled > baseQual[i].loc &&
        baseQual[i].locUncalled < baseQual[i+1].loc &&
        baseQual[i].flgUncalled )
    {
      valley[i] = tot_vals[baseQual[i].locUncalled];
    }
    else
    if( baseQual[i+1].locUncalled >= 0 &&
        baseQual[i+1].locUncalled > baseQual[i].loc &&
        baseQual[i+1].locUncalled < baseQual[i+1].loc &&
        baseQual[i+1].flgUncalled )
    {
      valley[i] = tot_vals[baseQual[i+1].locUncalled];
    }
  }

  win_wid = 7;
  win_edg = win_wid / 2;

  for( i = win_edg; i < numBase - win_edg; ++i )
  {
    lo = i - win_edg;
    hi = i + win_edg;

    notN = 1;
    min_peak = tot_vals[baseQual[lo].loc];
    max_valley = valley[lo];
    for( j = lo; j <= hi; ++j )
    {
      if( tot_vals[baseQual[j].loc] < min_peak )
      {
        min_peak = tot_vals[baseQual[j].loc];
      }
      if( j < hi && valley[j] > max_valley )
      {
        max_valley = valley[j];
      }
      if( baseQual[j].nuc == 4 )
      {
        notN = 0;
        break;
      }
    }

    baseQual[i].resPar7 = ( notN && min_peak > 0.0 ) ?
                          ( max_valley / min_peak ) :
                          10000.0;
  }

  for( i = 0; i < win_edg; ++i )
  {
    baseQual[i].resPar7           = baseQual[win_edg].resPar7;
    baseQual[numBase-i-1].resPar7 = baseQual[numBase-win_edg-1].resPar7;
  }

  /*
  ** Free memory.
  */
  free( valley );

  return( OK );
}


