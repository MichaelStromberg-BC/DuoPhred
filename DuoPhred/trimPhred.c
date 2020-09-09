/** trimPhred.c **/

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
** Trim the sequence using a modification of
** the method developed by Richard Mott at the
** Sanger Centre.
**
** We find the subsequence with the maximum score
** calculated as ( 0.05 - p_e ) of each base.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "phred.h"

#ifdef ANSI_C
int trimPhred( PhredData *phredData )
#else
int trimPhred( phredData )
PhredData *phredData;
#endif
{
  int i, n, tbg;
  int beg, end;
  int indices[128];
  int *qv;
  float probScore;
  float maxScore;
  float minPrbVal;
  Option *option;

  option = getOption();

  if( option->trimSetOption == 0 )
  {
    minPrbVal            = TRIM_MIN_PRB_VAL;
  }
  else
  {
    minPrbVal = option->trimSetValue;
  }

  phredData->trimSetValue = minPrbVal;

  n  = phredData->numBase[LCL];
  qv = phredData->baseQual[LCL];

  /*
  ** Find the maximum scoring subsequence.
  */
  beg       = 0;
  end       = 0;
  tbg       = 0;
  probScore = 0.0;
  maxScore  = 0.0;
  for( i = 0; i < n; ++i )
  {
    probScore += minPrbVal - (float)pow( (double)10.0, (double)( -qv[i] / 10.0 ) );
    if( probScore <= 0.0 )
    {
      probScore = 0.0;
      tbg = i + 1;
    }
    if( probScore > maxScore )
    {
      beg          = tbg;
      end          = i;
      maxScore     = probScore;
    }
  }

  /*
  ** Filter out very short sequences and sequences
  ** with low overall quality.
  */
  if( end - beg + 1 < TRIM_MIN_SEQ_LEN )
  {
    phredData->lftPhredTrim = -1;
    phredData->rhtPhredTrim = -1;
  }
  else
  {
    phredData->lftPhredTrim = beg;
    phredData->rhtPhredTrim = end;
  }

  /*
  ** Verbosity.
  */
  if( option->verboseOption && option->verboseLevel >= 16 )
  {
    fprintf( stderr,
             "trimPhred: trim bases with cutoff %4.2f: retain %d to %d (first base is zero)\n",
             minPrbVal,
             beg,
             end );
  }

  return( 0 );
}

