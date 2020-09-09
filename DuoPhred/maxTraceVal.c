/** maxVal.c **/


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

static FLOAT maxTraceVal_maxVal[4];

#ifdef ANSI_C
int setMaxVal( int tr_length, FLOAT **tr_vals )
#else
int setMaxVal( tr_length, tr_vals )
int tr_length;
FLOAT **tr_vals;
#endif
{
  int i;
  FLOAT max0, max1, max2, max3;
  FLOAT *trv0;
  FLOAT *trv1;
  FLOAT *trv2;
  FLOAT *trv3;

  trv0 = tr_vals[0];
  trv1 = tr_vals[1];
  trv2 = tr_vals[2];
  trv3 = tr_vals[3];
  max0 = tr_vals[0][0];
  max1 = tr_vals[1][0];
  max2 = tr_vals[2][0];
  max3 = tr_vals[3][0];
  for( i = 0; i < tr_length; ++i )
  {
    if( *trv0 > max0 ) max0 = *trv0;
    if( *trv1 > max1 ) max1 = *trv1;
    if( *trv2 > max2 ) max2 = *trv2;
    if( *trv3 > max3 ) max3 = *trv3;
    ++trv0;
    ++trv1;
    ++trv2;
    ++trv3;
  }
  maxTraceVal_maxVal[0] = max0;
  maxTraceVal_maxVal[1] = max1;
  maxTraceVal_maxVal[2] = max2;
  maxTraceVal_maxVal[3] = max3;

  return( 0 );
}
  

#ifdef ANSI_C
int getMaxVal( FLOAT *maxTraceVal_g )
#else
int getMaxVal( maxTraceVal_g )
FLOAT *maxTraceVal_g;
#endif
{
  maxTraceVal_g[0] = maxTraceVal_maxVal[0];
  maxTraceVal_g[1] = maxTraceVal_maxVal[1];
  maxTraceVal_g[2] = maxTraceVal_maxVal[2];
  maxTraceVal_g[3] = maxTraceVal_maxVal[3];

  return( 0 );
}



