/** report.c **/


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

#define MAXBIN		1024
#define BINWID		10


static int numBin;
static int bins[MAXBIN];

#ifdef ANSI_C
int initReport( void )
#else
int initReport()
#endif
{
  int i;

  for( i = 0; i < MAXBIN; ++i )
  {
    bins[i] = 0;
  }

  numBin = 0;

  return( 0 );
}



#ifdef ANSI_C
int makeReport( PhredData *phredData )
#else
int makeReport( phredData )
PhredData *phredData;
#endif
{
  int i;
  int sum;
  int bin;

  sum = 0;

  for( i = 0; i < phredData->numBase[LCL]; ++i )
  {
    if( phredData->baseQual[LCL][i] >= 20 )
    {
      ++sum;
    }
  }

  bin = (int)( (float)sum / (float)BINWID );

  if( bin < MAXBIN )
  {
    ++bins[bin];

    if( bin + 1 > numBin )
    {
      numBin = bin + 1;
    }
  }

  return( 0 );
}




#ifdef ANSI_C
int writeReport( Option *option )
#else
int writeReport( option )
Option *option;
#endif
{
  int i;
  FILE *fp;

  if( ( fp = fopen( option->qualReportFileName, "w+" ) ) == NULL )
  {
    return( ERROR );
  }

  fprintf( fp, "High Quality Bases\n" );
  fprintf( fp, "\n" );

  fprintf( fp, "number of high\n" );
  fprintf( fp, "quality bases     number of reads\n" );
  fprintf( fp, "-----             ---------------\n" );
  for( i = 0; i < numBin; ++i )
  {
    fprintf( fp, "%4d-%-4d          %-6d\n",
             i * BINWID, ( i + 1 ) * BINWID - 1, bins[i] );
  }
  fprintf( fp, "\n" );

  fclose( fp );

  return( 0 );

}
