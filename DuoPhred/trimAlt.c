/** trimAlt.c **/

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

#define MAX_ENZ_POS	100

#ifdef ANSI_C
int trimAlt( PhredData *phredData, char *enzName )
#else
int trimAlt( phredData, enzName )
PhredData *phredData;
char *enzName;
#endif
{
  int i, n, tbg;
  int beg, end;
  int enzSite;
  int indices[128];
  int *qv;
  float probScore;
  float qualScore;
  float maxQualScore;
  float maxScore;
  float minPrbVal;
  Option *option;

  option = getOption();

  if( option->trimSetOption == 0 )
  {
    minPrbVal = TRIM_MIN_PRB_VAL;
  }
  else
  {
    minPrbVal = option->trimSetValue;
  }

  phredData->trimSetValue = minPrbVal;

  n  = phredData->numBase[LCL];
  qv = phredData->baseQual[LCL];

  /*
  ** Search for enzyme sequence if requested.
  */
  enzSite = -1;
  if( strlen( enzName ) != 0 )
  {
    int nn;
    for( i = 0; i < 3; ++i )
    {
      nn = stringMatch( enzName, strlen( enzName ),
                        phredData->base[LCL], n, i, indices );
      if( nn > 0 )
      {
        if( indices[0] < MAX_ENZ_POS )
        {
          enzSite = indices[0] + strlen( enzName );
          break;
        }
      }
    }
  }

  /*
  ** Find the maximum scoring subsequence.
  */
  beg       = 0;
  end       = 0;
  tbg       = 0;
  probScore = 0.0;
  maxScore  = 0.0;
  qualScore = 0.0;
  maxQualScore = 0.0;
  for( i = 0; i < n; ++i )
  {
    probScore += minPrbVal - (float)pow( (double)10.0, (double)( -qv[i] / 10.0 ) );
    qualScore += (float)qv[i];
    if( probScore <= 0.0 )
    {
      probScore = 0.0;
      qualScore = 0.0;
      tbg = i + 1;
    }
    if( probScore > maxScore )
    {
      beg          = tbg;
      end          = i;
      maxScore     = probScore;
      maxQualScore = qualScore;
    }
  }

  /*
  ** Enzyme site is start if requested.
  */
  if( enzSite >= 0 && enzSite > beg )
  {
    beg = enzSite;
    if( beg > end )
    {
      beg = end;
    }
  }

  /*
  ** A correction required for consistency with ted and
  ** original phred trimming code.
  */
  ++end;
  if( end >= n )
  {
    end = n - 1;
  }

  /*
  ** Filter out very short sequences.
  */
  if( end - beg + 1 < TRIM_MIN_SEQ_LEN )
  {
    beg = 0;
    end = 0;
  }

  /*
  ** base[beg] is the first retained base and
  ** base[end] is the last retained base.  Store
  ** these values in a manner consistent with
  ** the `ted' display.
  */
  phredData->leftTrimPoint = beg;
  phredData->rghtTrimPoint = n - end;

  /*
  ** Verbosity.
  */
  if( option->verboseOption && option->verboseLevel >= 16 )
  {
    fprintf( stderr,
             "trimAlt: trim bases with cutoff %4.2f: retain %d to %d (first base is zero)\n",
             minPrbVal,
             beg,
             end - 1 );
  }

  return( 0 );
}

