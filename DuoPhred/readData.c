/** readData.c **/

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
*******************************************************************************
**                                                                           **
**    * readData.c                                                     *     **
**    * benefits from ideas in code written by LaDeana Hillier and     *     **
**    * Tim Gleeson.                                                   *     **
**                                                                           **
*******************************************************************************
*/

/*
** Note: try to have phred write a 'empty' output files
**       on error.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "phred.h"
#include "readType.h"
#include "readSCF.h"
#include "readABI.h"
#include "readESD.h"

/*
** Read sequence data.
*/
#ifdef ANSI_C
ChromatData *readData( char *filename, int *status )
#else
ChromatData *readData( filename, status )
char *filename;
int *status;
#endif
{
  int inType;
  int lstat;
  ChromatData *chromatData;
  Option *option;

  option = getOption();

  inType = readType( filename );
  if( inType == -1 )
  {
    fprintf( stderr, "error: file not found > %s\n", filename );
    *status = 1;
    return( NULL );
  }

  /*
  ** read data
  */ 
  switch( inType )
  {
    case SCFFormat:
      /*
      ** SCF
      */

      if( option->verboseOption && option->verboseLevel >= 16 )
      {
        fprintf( stderr,
                 "readData: SCF format file: %s\n",
                 filename );
      }

      chromatData = readSCF( filename, &lstat );
      *status = lstat;
      break;

    case ABIFormat:
      /*
      ** ABI
      */

      if( option->verboseOption && option->verboseLevel >= 16 )
      {
        fprintf( stderr,
                 "readData: ABI format file: %s\n",
                 filename );
      }

      chromatData = readABI( filename, &lstat );
      *status = lstat;
      break;

    case MD1Format:
      /*
      ** ESD format one (old ESD file)
      */

      if( option->verboseOption && option->verboseLevel >= 16 )
      {
        fprintf( stderr,
                 "readData: ESD format file: without extended header: %s\n",
                 filename );
      }

      chromatData = readESD( filename, MD1Format, &lstat );
      *status = lstat;
      break;

    case MD2Format:
      /*
      ** ESD format two (new ESD file)
      */

      if( option->verboseOption && option->verboseLevel >= 16 )
      {
        fprintf( stderr,
                 "readData: ESD format file: with extended header: %s\n",
                 filename );
      }

      chromatData = readESD( filename, MD2Format, &lstat );
      *status = lstat;
      break;

    default:
      fprintf( stderr, "  readData: unknown file type\n" );
      *status = 1;
      return( NULL );
  }

  return( chromatData );
}

