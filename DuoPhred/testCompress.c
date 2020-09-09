/** testCompress.c **/

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
#include <sys/types.h>
#include <sys/stat.h>
#include "phred.h"

/*
** Examine the first twp bytes of
** the file to determine whether
** it's compressed by gzip or
** compress.
*/
#ifdef ANSI_C
static int testCompressType( char *ifnm )
#else
static int testCompressType( ifnm )
char *ifnm;
#endif
{
  int type;
  unsigned char buf[3];
  FILE *fp;

  /*
  ** Read first three bytes of file.
  */
  type = -1;
  if( ( fp = fopen( ifnm, "rb" ) ) == NULL )
  {
    return( type );
  }

  if( fread( buf, 3, 1, fp ) != 1 )
  {
    fclose( fp );
    return( type );
  }

  fclose( fp );

  /*
  ** Are they magic?
  */
  if( ( buf[0] & 0xff ) == 0x1f &&
      ( buf[1] & 0xff ) == 0x8b )
  {
    type = GZ_COMPRESS;
  }
  else
  if( ( buf[0] == 0x1f && buf[1] == 0x9d ) ||
      ( buf[0] == 0x9d && buf[1] == 0x1f ) )
  {
    type = Z_COMPRESS;
  }
  else
  if( buf[0] == 'B' &&
      buf[1] == 'Z' &&
      buf[2] == 'h' )
  {
    type = BZ_COMPRESS;
  }
  else
  {
    type = NO_COMPRESS;
  }

  return( type );
}

/*
** Look for file as named as <ifnm> or
** <ifnm>.gz, <ifnm.Z>, or <ifnm>.bz2.
** Return type of compression (or none)
** and copy the name of the file found
** into ofnm.
*/
#ifdef ANSI_C
int testCompress( char *ifnm, char *ofnm )
#else
int testCompress( ifnm, ofnm )
char *ifnm;
char *ofnm;
#endif
{
  int len;
  int compress;
  int status;
  char string[2*PHRED_PATH_MAX];
  struct stat fs;

  compress = NO_COMPRESS;

  /*
  ** Does the filename have a compression suffix?
  */
  len = strlen( ifnm );
  if( ifnm[len-1] == 'Z' &&
      ifnm[len-2] == '.' )
  {
    compress = Z_COMPRESS;
    pstrcpy( ofnm, ifnm, PHRED_PATH_MAX );
    if( testCompressType( ofnm ) == compress )
    {
      return( compress );
    }
    return( NO_COMPRESS );
  }
  else
  if( ifnm[len-1] == 'z' &&
      ifnm[len-2] == 'g' &&
      ifnm[len-3] == '.' )
  {
    compress = GZ_COMPRESS;
    pstrcpy( ofnm, ifnm, PHRED_PATH_MAX );
    if( testCompressType( ofnm ) == compress )
    {
      return( compress );
    }
    return( NO_COMPRESS );
  }
  else
  if( ifnm[len-1] == '2' &&
      ifnm[len-2] == 'z' &&
      ifnm[len-3] == 'b' )
  {
    compress = BZ_COMPRESS;
    pstrcpy( ofnm, ifnm, PHRED_PATH_MAX );
    if( testCompressType( ofnm ) == compress )
    {
      return( compress );
    }
    return( NO_COMPRESS );
  }
  else
  {
    sprintf( string,
             "%.*s.gz",
             PHRED_PATH_MAX - 4,
             ifnm );
    status = stat( string, &fs );
    if( status == 0 )
    {
      compress = GZ_COMPRESS;
      pstrcpy( ofnm, string, PHRED_PATH_MAX );
      if( testCompressType( ofnm ) == compress )
      {
        return( compress );
      }
      return( NO_COMPRESS );
    }

    sprintf( string,
             "%.*s.Z",
             PHRED_PATH_MAX - 3,
             ifnm );
    status = stat( string, &fs );
    if( status == 0 )
    {
      compress = Z_COMPRESS;
      pstrcpy( ofnm, string, PHRED_PATH_MAX );
      if( testCompressType( ofnm ) == compress )
      {
        return( compress );
      }
      return( NO_COMPRESS );
    }

    sprintf( string,
             "%.*s.bz2",
             PHRED_PATH_MAX - 5,
             ifnm );
    status = stat( string, &fs );
    if( status == 0 )
    {
      compress = BZ_COMPRESS;
      pstrcpy( ofnm, string, PHRED_PATH_MAX );
      if( testCompressType( ofnm ) == compress )
      {
        return( compress );
      }
      return( NO_COMPRESS );
    }

    /*
    ** Is it a compressed file without a compression
    ** suffix?
    */
    compress = testCompressType( ifnm );
    if( compress == GZ_COMPRESS )
    {
      pstrcpy( ofnm, ifnm, PHRED_PATH_MAX );
      return( compress );
    }
    else
    if( compress == BZ_COMPRESS )
    {
      pstrcpy( ofnm, ifnm, PHRED_PATH_MAX );
      return( compress );
    }
    else
    if( compress == Z_COMPRESS )
    {
      fprintf( stderr,
               "testCompress: error: unable to uncompress file without .Z suffix\n" );
      pstrcpy( ofnm, ifnm, PHRED_PATH_MAX );
      return( NO_COMPRESS );
    }
  }
  pstrcpy( ofnm, ifnm, PHRED_PATH_MAX );

  return( NO_COMPRESS );
}

