/** readType.c **/

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
**    * readType.c                                                     *     **
**    * benefits from ideas in code written by LaDeana Hillier and     *     **
**    * Tim Gleeson.                                                   *     **
**                                                                           **
*******************************************************************************
*/

/* ---- Imports ---- */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "phred.h"

typedef struct
{
  int type;
  int offset;
  char *string;
  int len;
} Magic;

#define NUMMAG	4
Magic magic[NUMMAG] =
{
{ SCFFormat , 0,   ".scf",             4  } ,
{ ABIFormat , 0,   "ABIF",             4  } ,
{ ABIFormat , 128, "ABIF",             4  } ,
{ SCFFormat , 0,   "\234\330\300\000", 16 }
};

#ifdef ANSI_C
int readType( char *fn )
#else
int readType( fn )
char *fn;
#endif
{
  int i;
  int len;
  char buf[512];
  FILE *fp;

  /*
  ** Try to open file.
  */
  if( ( fp = fopen( fn, "rb" ) ) == NULL )
  {
    return( EEKFormat );
  }

  /*
  ** Check magic strings.
  */
  for( i = 0 ; i < NUMMAG; i++ )
  {
    if( fseek( fp, magic[i].offset, 0 ) == 0 )
    {
      len = magic[i].len;
      if( fread( buf, len, 1, fp ) == 1 )
      {
        if( strncmp( buf, magic[i].string, len ) == 0 )
        {
          fclose( fp );
          return( magic[i].type );
        }
      }
    }
  }

  fclose( fp );

  /*
  ** Check for Molecular Dynamics ESD file by calculating
  ** the number of data points based on the file/header
  ** sizes and comparing that to the value in the trailer,
  ** assuming we can read a trailer successfully.
  */
  i = isESD( fn );
  if( i == 1 )
  {
    return( MD1Format );
  }
  else
  if( i == 2 )
  {
    return( MD2Format );
  }

  return( UNKFormat );
}

