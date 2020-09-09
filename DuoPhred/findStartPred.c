/** findStartPred.c **/

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
** Locate the trace region where the peaks are evenly spaced.
** This is a good point to begin peak prediction.  Also
** set tr_qual array to indicate regions of the trace where the
** data quality degrades, in order that the predicted period 
** be held constant.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "phred.h"

#ifdef ANSI_C
int findStartPred( int tr_length, int npk, LocPeak *locPeak, int block_size, int *tr_qual )
#else
int findStartPred( tr_length, npk, locPeak, block_size, tr_qual )
int tr_length;
int npk;
LocPeak *locPeak;
int block_size;
int *tr_qual;
#endif
{
  int i, j, k;
  int jlo, jhi;
  int loclo;
  int startLoc;
  int initLoc;
  int iqual;
  FLOAT minn;
  Option *option;

  option = getOption();

  /*
  ** Locate highest quality trace point, avoiding points within
  ** block_size of the trace ends.
  */
  startLoc = 0;
  for( j = 0; j < npk; ++j )
  {
    if( locPeak[j].quality != -1 )
    {
      startLoc = locPeak[j].point;
      break;
    }
  }

  /*
  ** Few peaks?
  */
  if( npk < 20 )
  {
    for( i = 0; i < tr_length; ++i )
    {
      tr_qual[i] = -1;
    }
    initLoc = 0.25 * tr_length;
    return( initLoc );
  }

  /*
  ** Avoid ends.
  */
  jlo = 5;
  jhi = npk - 6;

  loclo = block_size > startLoc ? block_size : startLoc;
  while( locPeak[jlo].point < loclo &&
         locPeak[jlo].point < tr_length - block_size &&
         jlo < npk - 1 )
  {
    ++jlo;
  }
  while( locPeak[jhi].point >= tr_length - block_size && jhi > jlo )
  {
    --jhi;
  }
  minn = locPeak[jlo].stdev;
  k = jlo;
  for( j = jlo; j < jhi; ++j )
  {
    if( locPeak[j].quality == 0 &&
        locPeak[j].stdev < minn )
    {
      minn = locPeak[j].stdev;
      k = j;
    }
  }
  initLoc = locPeak[k].point;

  /*
  ** Assign trace quality used to control
  ** period variation during peak prediction.
  */
  jlo = 0;
  iqual = -1;
  for( i = 0; i < npk; ++i )
  {
    jhi = locPeak[i].point;
    for( j = jlo; j < jhi; ++j )
    {
      tr_qual[j] = iqual;
    }
    jlo = jhi;
    iqual = locPeak[i].quality;
  }
  if( jlo < tr_length - 1 )
  {
    for( j = jlo; j < tr_length; ++j )
    {
      tr_qual[j] = iqual;
    }
  }


  return( initLoc );
}

