/** polyData.c **/

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
** Functions related to writing data files for
** polymorphism study.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "phred.h"

#ifdef ANSI_C
int loadPolyData( Peak *first_peak, FLOAT *scale_traces, FLOAT **tr_vals, PolyData *polyData )
#else
int loadPolyData( first_peak, scale_traces, tr_vals, polyData )
Peak *first_peak;
FLOAT *scale_traces;
FLOAT **tr_vals;
PolyData *polyData;
#endif
{
  int i;
  int type, shift;
  int numBase;
  Peak *peak;
  Observed_peak *obs_peak;

  /*
  ** Count bases.
  */
  numBase = 0;
  for( peak = first_peak; peak; peak = peak->next )
  {
    ++numBase;
  }

  /*
  ** Initialize poly data.
  */
  polyData->called_nuc             = NULL;
  polyData->called_loc             = NULL;
  polyData->called_area            = NULL;
  polyData->called_relative_area   = NULL;
  polyData->uncalled_nuc           = NULL;
  polyData->uncalled_loc           = NULL;
  polyData->uncalled_area          = NULL;
  polyData->uncalled_relative_area = NULL;
  polyData->A_tr_val               = NULL;
  polyData->C_tr_val               = NULL;
  polyData->G_tr_val               = NULL;
  polyData->T_tr_val               = NULL;

  /*
  ** Allocate arrays.
  */
  polyData->alloc = 1;
  polyData->called_nuc             = (int *)ourMalloc( numBase * sizeof( int ) );
  polyData->called_loc             = (int *)ourMalloc( numBase * sizeof( int ) );
  polyData->called_area            = (FLOAT *)ourMalloc( numBase * sizeof( FLOAT ) );
  polyData->called_relative_area   = (FLOAT *)ourMalloc( numBase * sizeof( FLOAT ) );
  polyData->uncalled_nuc           = (int *)ourMalloc( numBase * sizeof( int ) );
  polyData->uncalled_loc           = (int *)ourMalloc( numBase * sizeof( int ) );
  polyData->uncalled_area          = (FLOAT *)ourMalloc( numBase * sizeof( FLOAT ) );
  polyData->uncalled_relative_area = (FLOAT *)ourMalloc( numBase * sizeof( FLOAT ) );
  polyData->A_tr_val               = (FLOAT *)ourMalloc( numBase * sizeof( FLOAT ) );
  polyData->C_tr_val               = (FLOAT *)ourMalloc( numBase * sizeof( FLOAT ) );
  polyData->G_tr_val               = (FLOAT *)ourMalloc( numBase * sizeof( FLOAT ) );
  polyData->T_tr_val               = (FLOAT *)ourMalloc( numBase * sizeof( FLOAT ) );

  /*
  ** Load scalars.
  */
  polyData->numBase = numBase;
  polyData->A_scale_factor = scale_traces[0];
  polyData->C_scale_factor = scale_traces[1];
  polyData->G_scale_factor = scale_traces[2];
  polyData->T_scale_factor = scale_traces[3];

  /*
  ** Load arrays.
  */
  i = 0;
  for( peak = first_peak; peak; peak = peak->next )
  {
    obs_peak = peak->obs_peak;
    if( obs_peak )
    {
      type = obs_peak->type;
      for( shift = 1; shift <= type; ++shift )
      {
        polyData->called_nuc[i]           = obs_peak->nuc;
        polyData->called_loc[i]           = obs_peak->split[type][shift];
        polyData->called_area[i]          = obs_peak->area;
        polyData->called_relative_area[i] = obs_peak->relative_area;
        polyData->A_tr_val[i]             = tr_vals[0][polyData->called_loc[i]];
        polyData->C_tr_val[i]             = tr_vals[1][polyData->called_loc[i]];
        polyData->G_tr_val[i]             = tr_vals[2][polyData->called_loc[i]];
        polyData->T_tr_val[i]             = tr_vals[3][polyData->called_loc[i]];
        if( peak->best_uncalled_peak )
        {
          polyData->uncalled_nuc[i]           = peak->best_uncalled_peak->nuc;
          polyData->uncalled_loc[i]           = peak->best_uncalled_peak->split[1][1];
          polyData->uncalled_area[i]          = peak->best_uncalled_peak->area;
          polyData->uncalled_relative_area[i] = peak->best_uncalled_peak->relative_area;
        }
        else
        {
          polyData->uncalled_nuc[i]           = 4;
          polyData->uncalled_loc[i]           = -1;
          polyData->uncalled_area[i]          = -1.0;
          polyData->uncalled_relative_area[i] = -1.0;
        }
        if( shift != type )
        {
          peak = peak->next;
        }
        ++i;
      }
    }
    else
    {
      polyData->called_nuc[i]           = 4;
      polyData->called_loc[i]           = (int)peak->pred_location;
      polyData->called_area[i]          = -1.0;
      polyData->called_relative_area[i] = -1.0;
      polyData->A_tr_val[i]             = tr_vals[0][polyData->called_loc[i]];
      polyData->C_tr_val[i]             = tr_vals[1][polyData->called_loc[i]];
      polyData->G_tr_val[i]             = tr_vals[2][polyData->called_loc[i]];
      polyData->T_tr_val[i]             = tr_vals[3][polyData->called_loc[i]];
      if( peak->best_uncalled_peak )
      {
        polyData->uncalled_nuc[i]           = peak->best_uncalled_peak->nuc;
        polyData->uncalled_loc[i]           = peak->best_uncalled_peak->split[1][1];
        polyData->uncalled_area[i]          = peak->best_uncalled_peak->area;
        polyData->uncalled_relative_area[i] = peak->best_uncalled_peak->relative_area;
      }
      else
      {
        polyData->uncalled_nuc[i]           = 4;
        polyData->uncalled_loc[i]           = -1;
        polyData->uncalled_area[i]          = -1.0;
        polyData->uncalled_relative_area[i] = -1.0;
      }
      ++i;
    }
  }

  return( OK );
}


#ifdef ANSI_C
int freePolyData( PolyData *polyData )
#else
int freePolyData( polyData )
PolyData *polyData;
#endif
{
  if( polyData->called_nuc != NULL )
  {
    ourFree( (char *)polyData->called_nuc );
  }
  if( polyData->called_loc != NULL )
  {
    ourFree( (char *)polyData->called_loc );
  }
  if( polyData->called_area != NULL )
  {
    ourFree( (char *)polyData->called_area );
  }
  if( polyData->called_relative_area != NULL )
  {
    ourFree( (char *)polyData->called_relative_area );
  }
  if( polyData->uncalled_nuc != NULL )
  {
    ourFree( (char *)polyData->uncalled_nuc );
  }
  if( polyData->uncalled_loc != NULL )
  {
    ourFree( (char *)polyData->uncalled_loc );
  }
  if( polyData->uncalled_area != NULL )
  {
    ourFree( (char *)polyData->uncalled_area );
  }
  if( polyData->uncalled_relative_area != NULL )
  {
    ourFree( (char *)polyData->uncalled_relative_area );
  }
  if( polyData->A_tr_val != NULL )
  {
    ourFree( (char *)polyData->A_tr_val );
  }
  if( polyData->C_tr_val != NULL )
  {
    ourFree( (char *)polyData->C_tr_val );
  }
  if( polyData->G_tr_val != NULL )
  {
    ourFree( (char *)polyData->G_tr_val );
  }
  if( polyData->T_tr_val != NULL )
  {
    ourFree( (char *)polyData->T_tr_val );
  }

  return( OK );
}



#ifdef ANSI_C
int writePolyData( char *filename, char *seqname, int numPoint, PolyData *polyData )
#else
int writePolyData( filename, seqname, numPoint, polyData )
char *filename;
char *seqname;
int numPoint;
PolyData *polyData;
#endif
{
  int i;
  FLOAT minScale;
  FLOAT scale_traces[4];
  FILE *fp;

  /*
  ** Open output file.
  */
  if( ( fp = fopen( filename, "w+" ) ) == NULL )
  {
    fprintf( stderr, "writePolyData: error: unable to open file %s\n", filename );
    return( ERROR );
  }

  if( polyData->called_nuc != NULL && numPoint > 0 )
  {
    /*
    ** Find minimum scale factor.
    */
    scale_traces[0] = polyData->A_scale_factor;
    scale_traces[1] = polyData->C_scale_factor;
    scale_traces[2] = polyData->G_scale_factor;
    scale_traces[3] = polyData->T_scale_factor;
    for( i = 0; i < 4; ++i )
    {
      if( scale_traces[i] > 0.0 )
      {
        minScale = scale_traces[i];
        break;
      }
    }
    for( i = 0; i < 4; ++i )
    {
      if( scale_traces[i] > 0.0 && scale_traces[i] < minScale )
      {
        minScale = scale_traces[i];
      }
    }

    /*
    ** Write header.
    */
    fprintf( fp, "%s %f %f %f %f %f\n",
                 seqname,
                 minScale,
                 polyData->A_scale_factor,
                 polyData->C_scale_factor,
                 polyData->G_scale_factor,
                 polyData->T_scale_factor );

    /*
    ** Write data.
    */
    for( i = 0; i < polyData->numBase; ++i )
    {
      fprintf( fp, "%c  %d  %f  %f  %c  %d  %f  %f  %f  %f  %f  %f",
               "ACGTN"[polyData->called_nuc[i]],
               polyData->called_loc[i],
               polyData->called_area[i],
               polyData->called_relative_area[i],
               "ACGTN"[polyData->uncalled_nuc[i]],
               polyData->uncalled_loc[i],
               polyData->uncalled_area[i],
               polyData->uncalled_relative_area[i],
               polyData->A_tr_val[i],
               polyData->C_tr_val[i],
               polyData->G_tr_val[i],
               polyData->T_tr_val[i] );
      if( fwrite( "\n", 1, 1, fp ) != 1 )
      {
        fprintf( stderr,
                 "writePolyData: error: unable to write to %s\n",
                 filename );
        fclose( fp );
        delFile( filename );
        return( ERROR );
      }
    }

  }

  /*
  ** Close file.
  */
  fclose( fp );

  return( OK );
}

