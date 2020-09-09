/** makePeakSignal.c **/

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
** Generate "square" waveform with a period that corresponds
** as closely as possible to the period of the observed peaks.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "phred.h"

#ifdef ANSI_C
int makePeakSignal( int tr_length, int npk, LocPeak *locPeak, FLOAT *waveForm )
#else
int makePeakSignal( tr_length, npk, locPeak, waveForm )
int tr_length;
int npk;
LocPeak *locPeak;
FLOAT *waveForm;
#endif
{
  int i, j;
  int loc;
  int jlo, jhi;
  Option *option;

  option = getOption();

  /*
  ** Create a "square" wave for use by the
  ** make_predicted_peaks() module.
  */
  if( waveForm != NULL )
  {
    for( i = 0; i < tr_length; ++i )
    {
      waveForm[i] = -1.0;
    }
    /*
    ** first peak...
    */
    jlo = locPeak[0].point -
          ( locPeak[1].point - locPeak[0].point ) * 0.25;
    jhi = locPeak[0].point +
          ( locPeak[1].point - locPeak[0].point ) * 0.25;
    if( jlo < 0 )
    {
      jlo = 0;
    }
    for( j = jlo + 1; j <= jhi; ++j )
    {
      waveForm[j] = 1.0;
    }
    /*
    ** all but first and last peaks...
    */
    for( i = 1; i < npk - 1; ++i )
    {
      loc = locPeak[i].point;
      jlo = loc - ( loc - locPeak[i-1].point ) * 0.25;
      jhi = loc + ( locPeak[i+1].point - loc ) * 0.25;
      for( j = jlo + 1; j <= jhi; ++j )
      {
        waveForm[j] = 1.0;
      }
    }
    /*
    ** last peak...
    */
    jlo = locPeak[npk-1].point -
          ( locPeak[npk-1].point - locPeak[npk-2].point ) * 0.25;
    jhi = locPeak[npk-1].point +
          ( locPeak[npk-1].point - locPeak[npk-2].point ) * 0.25;
    if( jhi >= tr_length )
    {
      jhi = tr_length - 1;
    }
    for( j = jlo + 1; j <= jhi; ++j )
    {
      waveForm[j] = 1.0;
    }
  }


  return( OK );
}
