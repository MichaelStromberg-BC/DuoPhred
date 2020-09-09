/** evaluateTrace.c **/

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
** Evaluate the trace quality using several measures.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "phred.h"

/*
** cmpLocPeak
**
** Used by qsort() to order peaks by position in trace.
*/
#ifdef ANSI_C
static int cmpLocPeak( LocPeak *peak1, LocPeak *peak2 )
#else
static int cmpLocPeak( peak1, peak2 )
LocPeak *peak1;
LocPeak *peak2;
#endif
{
  if( peak1->point == peak2->point )
  {
    return( 0 );
  }
  else
  if( peak1->point < peak2->point )
  {
    return( -1 );
  }
  else
  {
    return( 1 );
  }
}



/*
** findCenter
**
** Find center of peak between loPnt and hiPnt.
**
*/
#ifdef ANSI_C
static int findCenter( int loPnt, int hiPnt, FLOAT *tr )
#else
static int findCenter( loPnt, hiPnt, tr )
int loPnt;
int hiPnt;
FLOAT *tr;
#endif
{
  int i;
  int loc;
  int iflag;
  register FLOAT ndrv;
  register FLOAT odrv;

  /*
  ** Look for zero crossing in first derivative.
  */
  iflag = 0;
  loc = loPnt;
  odrv = tr[loPnt+1] - tr[loPnt];
  for( i = loPnt + 1; i <= hiPnt; ++i )
  {
    ndrv = tr[i+1] - tr[i];
    if( odrv > 0.0 && ndrv <= 0.0 )
    {
      loc = i;
      iflag = 1;
      break;
    }
    odrv = ndrv;
  }

  /*
  ** Assume a shoulder in absence of a maximum
  ** and choose a midpoint between inflection
  ** points. (Searching for a minimum in second
  ** derivative didn't do better.)
  */
  if( iflag == 0 )
  {
    loc = - ( loPnt + hiPnt ) * 0.5;
  }

  return( loc );
}


/*
** findSubsumed
**
** Purpose:
**
** Look for truncated peaks that probably consist
** of several peaks and split them up.
*/
#ifdef ANSI_C
static int findSubsumed( int tr_length, FLOAT **trv, FLOAT *max, LocPeak *locPeak, int *npk, int firstTrunc, int lastTrunc )
#else
static int findSubsumed( tr_length, trv, max, locPeak, npk, firstTrunc, lastTrunc )
int tr_length;
FLOAT **trv;
FLOAT *max;
LocPeak *locPeak;
int *npk;
int firstTrunc;
int lastTrunc;
#endif
{
  int i, j, n;
  int ilo, ihi;
  int loclo;
  int lochi;
  int iflag;
  int base;
  FLOAT var;
  FLOAT sum;
  FLOAT mean;
  FLOAT stddev;
  FLOAT width;

  ilo = 0;
  while( locPeak[ilo].point < firstTrunc &&
         ilo < ( *npk - 1 ) )
  {
    ++ilo;
  }

  ihi = ilo;
  while( locPeak[ihi].point < lastTrunc &&
         ihi < ( *npk - 1 ) )
  {
    ++ihi;
  }

  if( ilo == 0 )
  {
    /*
    ** Handle this case separately.
    ** (Ignore for now.)
    */
    ++ilo;
  }

  if( ihi == *npk - 1 )
  {
    /*
    ** Handle this case separately.
    ** (Ignore for now.)
    */
    --ihi;
  }

  for( i = ilo; i <= ihi; ++i )
  {
    /*
    ** Is this peak truncated?
    */
    if( trv[locPeak[i].base][locPeak[i].point] == max[locPeak[i].base] )
    {
      /*
      ** Find adjacent peak locations.
      */
      loclo = locPeak[i-1].point;
      lochi = locPeak[i+1].point;

      /*
      ** Search for nearest downstream set of four "good" peaks
      ** with which to infer the period.  I imagine this is not
      ** especially critical because I want to fill in peaks
      ** only where especially large gaps exists.
      ** 
      */
      iflag = 0;
      for( j = i; j < *npk - 4; ++j )
      {
        /*
        ** Do not wander too far away.
        */
        if( j > i + 30 )
        {
          iflag = -1;
          break;
        } 
        sum = ( locPeak[j+1].point - locPeak[j].point   ) +
              ( locPeak[j+2].point - locPeak[j+1].point ) +
              ( locPeak[j+3].point - locPeak[j+2].point ) +
              ( locPeak[j+4].point - locPeak[j+3].point );
        mean = sum / 4.0;

        /*
        ** Estimate variation amongst four peaks.
        */
        var = (FLOAT)( locPeak[j+1].point - locPeak[j].point   ) - mean;
        sum = var * var;
        var = (FLOAT)( locPeak[j+2].point - locPeak[j+1].point ) - mean;
        sum += var * var;
        var = (FLOAT)( locPeak[j+3].point - locPeak[j+2].point ) - mean;
        sum += var * var;
        var = (FLOAT)( locPeak[j+4].point - locPeak[j+3].point ) - mean;
        sum += var * var;

        stddev = (FLOAT)sqrt( sum / 3.0 );

        if( stddev / mean < 0.10 )
        {
          iflag = 1;
          break;
        }
      }
      if( iflag == 1 )
      {
        /*
        ** If the gap between two adjacent peaks is almost three times the
        ** apparent period or more, there are probably subsumed peaks.
        */
        if( (FLOAT)( lochi - loclo ) / mean > 2.8 )
        {
          /*
          ** Estimate number of subsumed peaks.
          */
          n =  (int)( (FLOAT)( lochi - loclo ) / mean + 0.2 ) - 1;

          /*
          ** Make space in the locPeak array.
          */
          for( j = *npk - 1; j > i; --j )
          {
            memcpy( &(locPeak[j+n-1]), &(locPeak[j]), sizeof( LocPeak ) );
          }

          /*
          ** Estimate positions of subsumed peaks and store them.
          */
          base = locPeak[i].base;
          width = (FLOAT)( lochi - loclo ) / (FLOAT)( n + 1 );
          for( j = 0; j < n; ++j )
          {
            locPeak[i+j].point = locPeak[i+j-1].point + width + 0.5;
            locPeak[i+j].base  = base;
            locPeak[i+j].type  = 2;
          }

          *npk += n - 1;
          i    += n - 1;
          ihi  += n - 1;
        }
      }
    }
  }

  return( 0 );
}


/*
** findPeak
**
** Locate peaks in trace.
**
** Look for inflection points.  When the
** second derivative goes from positive
** to negative, record the point assuming
** a peak occurs ahead.  When the second
** derivative goes from negative to positive,
** try to find a maximum between the two
** inflection points.  Return the maximum
** if it exists; otherwise, return the
** midpoint between the inflection points.
** Also try to find the center of clipped
** peaks.
** 
*/
#ifdef ANSI_C
static LocPeak *findPeak( int tr_length, FLOAT **trv, FLOAT *max, int *npk )
#else
static LocPeak *findPeak( tr_length, trv, max, npk )
int tr_length;
FLOAT **trv;
FLOAT *max;
int *npk;
#endif
{
  register int i;
  int j, k;
  register int loc;
  int ilead;
  int loclo, lochi;
  int iskip;
  int firstTrunc;
  int lastTrunc;
  int shoulder;
  int width;
  int truncated;
  register FLOAT val;
  register FLOAT ndrv, odrv;
  register FLOAT *ptrv;
  LocPeak *locPeak;

  locPeak = NULL;
  locPeak = (LocPeak *)ourMalloc( tr_length * sizeof( LocPeak ) );

  /*
  ** Loop through traces.
  */
  *npk = 0;
  firstTrunc = tr_length;
  lastTrunc = 0;
  for( j = 0; j < 4; ++j )
  {
    /*
    ** Initialize old derivative and
    ** location of leading peak shoulder.
    ** (ilead also serves as a flag indicating
    ** that a leading shoulder was found.)
    */
    odrv = 0.0;
    ilead = 0;
    ptrv = trv[j];

    /*
    ** Run along trace.
    */
    for( i = 1; i < tr_length - 1; ++i )
    {
      /*
      ** Reduce array indexing a bit.
      */
      val = ptrv[i];
      /*
      ** Second derivative.
      */
      ndrv = ( ptrv[i+1] - val ) - ( val - ptrv[i-1] );
      if( odrv > 0.0 && ndrv <= 0.0 )
      {
        if( i < tr_length - 10 &&
            ( ( ptrv[i+2] - ptrv[i+1] ) - ( ptrv[i+1] - ptrv[i] ) ) > 0.0 )
        {
          continue;
        }
        /*
        ** This is an inflection point on a leading shoulder.
        */
        ilead = i;
      }
      else
      if( odrv < 0.0 && ndrv >= 0.0 )
      {
        if( i < tr_length - 10 &&
            ( ( ptrv[i+2] - ptrv[i+1] ) - ( ptrv[i+1] - ptrv[i] ) ) < 0.0 )
        {
          continue;
        }
        /*
        ** This is an inflection point on a trailing shoulder.
        */
        if( ilead )
        {
          /*
          ** Find peak maximum if it exists, if not, find
          ** midpoint between ilead and i.
          */
          loc = findCenter( ilead, i, ptrv );

          shoulder = ( loc < 0 );
          width = i - ilead;
          loc = abs( loc );

          /*
          ** Check for truncated peak.  If so, find midpoint
          ** in run of maximum value trace points (assumes
          ** symmetry).
          */
          if( ptrv[loc] == max[j] )
          {
            truncated = 1;
            loclo = loc;
            while( ptrv[loclo] == max[j] &&
                   loclo > 0 )
            {
              --loclo;
            }
            lochi = loc;
            while( ptrv[lochi] == max[j] &&
                   lochi < tr_length - 1 )
            {
              ++lochi;
            }
            loc = ( loclo + lochi ) * 0.5 + 0.5;
            /*
            ** How many peaks are truncated at this position?
            */
            for( k = 0; k < 4; ++k )
            {
              if( trv[k][loc] == max[k] && k != j )
              {
                ++truncated;
              }
            }
          }
          else
          {
            truncated = 0;
          }

          /*
          ** Record this peak if it is the highest
          ** amplitude of the four traces at this
          ** point...
          */
          iskip = 0;
          val = ptrv[loc];
          /*
          ** (This could be improved I imagine.)
          */
          for( k = 0; k < 4; ++k )
          {
            if( k != j &&
                val < trv[k][loc] )
            {
              iskip = 1;
              break;
            }
          }
          /*
          ** ... and it is at least 10% the height
          ** of the previous peak.
          */
          if( iskip == 0 && *npk > 0 )
          {
            if( !locPeak[*npk-1].truncated &&
                val < 0.10 * trv[locPeak[*npk-1].base][locPeak[*npk-1].point] )
            {
              iskip = 1;
            }
          }
          if( iskip == 0 )
          {
            /*
            ** Record peak.
            */
            locPeak[*npk].point = loc;
            locPeak[*npk].base = j;
            locPeak[*npk].type = 1;
            locPeak[*npk].shoulder = shoulder;
            locPeak[*npk].width = width;
            locPeak[*npk].truncated = truncated;
            /*
            ** Record positions of first and last
            ** truncated peaks.
            */
            if( truncated )
            {
              if( loc < firstTrunc )
              {
                firstTrunc = loc;
              }
              if( loc > lastTrunc && loc < tr_length - 1000 )
              {
                lastTrunc = loc;
              }
            }
            ++(*npk);
          }
          ilead = 0;
        }
      }
      odrv = ndrv;
    }
  }

  /*
  ** Sort peaks by location.
  */
#ifdef ANSI_C
  qsort( locPeak, (size_t) *npk, sizeof( LocPeak ),
         (int (*)(const void *, const void *))cmpLocPeak );
#else
  qsort( locPeak, (size_t) *npk, sizeof( LocPeak ), cmpLocPeak );
#endif

  /*
  ** Look for peaks subsumed under a truncated envelope.
  */
  if( firstTrunc != tr_length && lastTrunc != 0 )
  {
    findSubsumed( tr_length, trv, max, locPeak, npk, firstTrunc, lastTrunc );
  }

  return( locPeak );
}



/*
** evaluateTrace
**
** Evaluate traces.
**
** Locate peaks in traces.  Calculate period of peaks within
** block_size ranges of each peak, and the standard deviation
** of the period (normalized to the period). Calculate spcRatio
** in seven peak windows.  Locate start and end of "good" trace.
*/
#ifdef ANSI_C
LocPeak *evaluateTrace( int tr_length, FLOAT **tr_vals,
                        FLOAT *tot_vals, int block_size,
                        int *numLocPeak,
                        int *begGoodTrace, int *endGoodTrace )
#else
LocPeak *evaluateTrace( tr_length, tr_vals, tot_vals, block_size,
                        numLocPeak,
                        begGoodTrace, endGoodTrace )
int tr_length;
FLOAT **tr_vals;
FLOAT *tot_vals;
int block_size;
int *numLocPeak;
int *begGoodTrace;
int *endGoodTrace;
#endif
{
  int owid;
  int nwid;
  int otrunc;
  int ntrunc;
  register int i, j, k;
  register int jlo, jhi;
  int loc, loclo, lochi;
  int npk;
  int pnt;
  register FLOAT sum;
  FLOAT fdif;
  FLOAT mean;
  FLOAT var;
  FLOAT wgt;
  FLOAT div;
  FLOAT maxx;
  FLOAT minn;
  FLOAT max[4];
  char line[2*PHRED_PATH_MAX];
  register LocPeak *locPeak;
  Option *option;

  char *ourMalloc();

  option = getOption();

  /*
  ** Get the maximum trace values for tr_vals.
  */
  getMaxVal( max );

  /*
  ** locate peaks in each trace
  */
  locPeak = findPeak( tr_length, tr_vals, max, &npk );

  if( npk < 8 )
  {
    if( option->tagOption )
    {
      sprintf( line,
               "PROCESSING_ERROR: %.*s: unable to call bases: too few peaks\n",
               PHRED_PATH_MAX,
               option->inFileName[option->curInFile] );
    }
    else
    {
      sprintf( line,
               "unable to call bases: too few peaks\n" );

    }
    fprintf( stderr, "%s", line );
    if( option->log )
    {
      writeLog( line );
    }

    if( locPeak != NULL )
    {
      ourFree( (char *)locPeak );
    }
    return( NULL );
  }

  /*
  ** Calculate spacing ratio.
  */
  for( i = 0; i < ( npk - 1 ); ++i )
  {
    locPeak[i].space = locPeak[i+1].point - locPeak[i].point;
  }
  locPeak[npk-1].space = locPeak[npk-2].space;

  for( i = 3; i < ( npk - 2 ); ++i )
  {
    minn = locPeak[i-3].space;
    maxx = minn;
    for( j = -2; j <= 2; ++j )
    {
      if( locPeak[i+j].space < minn )
      {
        minn = locPeak[i+j].space;
      }
      if( locPeak[i+j].space > maxx )
      {
        maxx = locPeak[i+j].space;
      }
    }
    locPeak[i].spcRatio = ( minn > 0.0 ) ? maxx / minn : 100.0;
  }
  locPeak[0].spcRatio = locPeak[3].spcRatio;
  locPeak[1].spcRatio = locPeak[3].spcRatio;
  locPeak[2].spcRatio = locPeak[3].spcRatio;
  locPeak[npk-2].spcRatio = locPeak[npk-3].spcRatio;
  locPeak[npk-1].spcRatio = locPeak[npk-3].spcRatio;

  /*
  ** Attempt to identify at which point useful trace begins and
  ** at which point it ends.
  */
  *begGoodTrace = 0;
  *endGoodTrace = tr_length;

  /*
  ** First try to find where the useful trace begins by searching
  ** for the "primer peak", which I believe saturates the four
  ** traces at the same time and occurs in the first half
  ** of the trace.
  */
  for( i = 0; i < npk; ++i )
  {
    if( locPeak[i].truncated == 4 && locPeak[i].point < 0.5 * tr_length )
    {
      *begGoodTrace = locPeak[i].point;
    }
  }

  /*
  ** Try to find where the useful trace ends by comparing
  ** the widths of adjacent peaks.  I believe that the
  ** the peaks widths vary strongly from peak to peak
  ** when the peaks fall to amplitudes near the level of
  ** noise.
  */
  owid = locPeak[0].width;
  otrunc = locPeak[0].truncated;
  for( i = *begGoodTrace; i < npk; ++i )
  {
    if( !(locPeak[i].shoulder) )
    {
      nwid = locPeak[i].width;
      ntrunc = locPeak[i].truncated;
      fdif = (FLOAT)abs( nwid - owid ) / (FLOAT)( ( nwid < owid ) ? nwid : owid );
      if( fdif > 3.0 && !otrunc && !ntrunc )
      {
        *endGoodTrace = locPeak[i].point;
        break;
      }
      owid = nwid;
      otrunc = ntrunc;
  
    }
  }


  /*
  ** Calculate the standard deviation of the period
  ** in block_size ranges around each peak.  Record
  ** the standard deviation as a fraction of the
  ** mean period in that range of peaks.
  */

  /*
  ** initialize result array...
  */
  for( i = 0; i < npk; ++i )
  {
    locPeak[i].stdev = 1.0;
  }

  /*
  ** loop through peaks...
  */
  for( i = 1; i < npk - 1; ++i )
  {
    /*
    ** find peaks about block_size around the
    ** peak of interest...
    */
    loc = locPeak[i].point;
    jlo = i;
    while( ( loc - locPeak[jlo].point ) < ( block_size / 2 ) &&
           jlo > 1 )
    {
      --jlo;
    }
    loclo = locPeak[jlo].point;

    jhi = i;
    while( ( locPeak[jhi].point - loc ) < ( block_size / 2 ) &&
           jhi < ( npk - 2 ) )
    {
      ++jhi;
    }
    lochi = locPeak[jhi].point;

    /*
    ** calculate mean period...
    */
    k = 0;
    sum = 0.0;
    for( j = jlo; j <= jhi; ++j )
    {
      sum += (FLOAT)( locPeak[j+1].point - locPeak[j].point );
      ++k;
    }

    /*
    ** check whether there is a significant number of peaks in
    ** this block...
    */
    if( k > 3 )
    {
      mean = sum / (FLOAT)k;
    }
    else
    {
      locPeak[i].stdev = 1000.0;
      continue;
    }

    /*
    ** calculate standard deviation...
    */
    div = 0.0;
    sum = 0.0;
    for( j = jlo; j < i; ++j )
    {
      pnt = locPeak[j].point;
      wgt = ( loc - loclo > 0 ) ?
            (FLOAT)( pnt - loclo ) / (FLOAT)( loc - loclo ) : 1.0;
      var = mean - (FLOAT)( pnt - locPeak[j-1].point );
      sum += ( var * var ) * wgt;
      div += wgt;
    }
    for( j = i; j <= jhi; ++j )
    {
      pnt = locPeak[j].point;
      wgt = ( lochi - loc > 0 ) ?
            (FLOAT)( lochi - pnt ) / (FLOAT)( lochi - loc ) :
            1.0;
      var = mean - (FLOAT)( locPeak[j].point - locPeak[j-1].point );
      sum += ( var * var ) * wgt;
      div += wgt;
    }
    locPeak[i].stdev = (FLOAT)sqrt( (double)( sum / div ) ) / mean;
  }

  /*
  ** Identify "bad" peaks.  Identify peaks
  ** lying before the trace maxima and
  ** peaks lying in areas with large peak
  ** spacing variability.
  */

  /*
  ** chop off low signal trace at start...
  */
  for( i = 0; i < npk; ++i )
  {
    locPeak[i].quality = 0;
  }

/*
  for( i = 0; i < npk; ++i )
  {
    k = locPeak[i].point;
    j  = locPeak[i].base;
    if( tr_vals[j][k] >= 0.70 * max[j] )
    {
      break;
    }
    locPeak[i].quality = -1;
  }
*/

  /*
  ** block out bad trace based on
  ** period variability...
  */
  for( i = 0; i < npk; ++i )
  {
    if( locPeak[i].quality == 0 &&
        locPeak[i].stdev > 0.45 )
    {
      locPeak[i].quality = 1;
    }
  }

  *numLocPeak = npk;

  return( locPeak );
}


