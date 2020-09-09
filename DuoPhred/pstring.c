/** pstring.c **/

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
#include "pstring.h"

/*
** Note: these functions save up to dst_len characters including
**       the terminating null character, which will always be
**       present. Thus, if src is longer than dst_len, the
**       resulting dst will contain ( dstlen - 1 ) characters
**       of src and a terminating null.
*/

#ifdef ANSI_C
char *pstrcpy( char *dst, char *src, int dst_len )
#else
char *pstrcpy( dst, src, dst_len )
char *dst;
char *src;
int dst_len;
#endif
{
  char *optr;
  char *eob;

  eob  = dst + dst_len - 1;
  optr = dst;

  while( *src != '\0' )
  {
    *dst = *src;

    if( dst >= eob )
    {
      fprintf( stderr,
               "pstrcpy: source string length exceeds destination buffer size: truncating\n" );
      break;
    }

    ++dst;
    ++src;
  }

  *dst = '\0';

  return( optr );
}

#ifdef ANSI_C
char *pstrcat( char *dst, char *src, int dst_len )
#else
char *pstrcat( dst, src, dst_len )
char *dst;
char *src;
int dst_len;
#endif
{
  char *optr;
  char *eob;

  eob  = dst + dst_len - 1;
  optr = dst;

  while( *dst != '\0' )
  {
    ++dst;
  }

  if( dst >= eob )
  {
    fprintf( stderr,
             "pstrcat: destination string exceeds destination buffer length: skip concatenation\n" );
    return( optr );
  }

  while( *src != '\0' )
  {
    *dst = *src;

    if( dst >= eob )
    {
      fprintf( stderr,
               "pstrcat: source string length exceeds destination buffer size: truncating\n" );
      break;
    }

    ++dst;
    ++src;
  }

  *dst = '\0';

  return( optr );
}

