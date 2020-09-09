/** readSCFUtil.c **/

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
**    * readSCF.c                                                      *     **
**    * benefits from ideas in code written by LaDeana Hillier and     *     **
**    * Tim Gleeson.                                                   *     **
**                                                                           **
*******************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <sys/types.h>
#include "typeDef.h"
#include "rwUtil.h"
#include "chromatData.h"
#include "freeChromatData.h"
#include "readSCF.h"

#ifdef ANSI_C
CommentEntry *parseSCFComment( int len, char *comment, int *numEntry )
#else
CommentEntry *parseSCFComment( len, comment, numEntry )
int len;
char *comment;
int *numEntry;
#endif
{
  int i, n;
  int istate;
  int numByte;
  CommentEntry *ce;

  /*
  ** Count carriage returns (the first
  ** comment does not begin with a
  ** linefeed.
  */
  n = 1;
  for( i = 0; i < len; ++i )
  {
    if( comment[i] == '\n' )
    {
      ++n;
    }
  }
  ++n;

  /*
  ** Allocate and initialize memory.
  */
  numByte = n * sizeof( CommentEntry );
  ce = (CommentEntry *)malloc( numByte );
  if( ce == NULL )
  {
    fprintf( stderr, "parseSCFComment: error: unable to allocate memory\n" );
    return( NULL );
  }
  for( i = 0; i < n; ++i )
  {
    ce[i].label = NULL;
    ce[i].value = NULL;
    ce[i].used  = 0;
  }

  n = 0;
  istate = 0;
  for( i = 0; i < len; ++i )
  {
    if( istate == 0 )
    {
      /*
      ** Wait for start of label.
      */
      if( comment[i] != ' ' &&
          comment[i] != '\n' )
      {
        ce[n].label = &(comment[i]);
        istate = 1;
      }
    }
    else
    if( istate == 1 )
    {
      /*
      ** Wait for end of label.
      */
      if( comment[i] == ' ' ||
          comment[i] == '=' )
      {
        comment[i] = '\0';
        istate = 2;
      }
      else
      if( comment[i] == '\n' ||
          comment[i] == '\0' )
      {
        comment[i] = '\0';
        istate = 0;
      }
    }
    else
    if( istate == 2 )
    {
      /*
      ** Wait for start of value.
      */
      if( comment[i] == '\n' )
      {
        /*
        ** Discard NULL comments.
        */
        comment[i] = '\0';
        istate = 0;
      }
      else
      if( comment[i] != ' ' &&
          comment[i] != '=' )
      {
        ce[n].value = &(comment[i]);
        istate = 3;
      }
      else
      if( comment[i] == '\0' )
      {
        comment[i] = '\0';
        istate = 0;
      }
    }
    else
    if( istate == 3 )
    {
      /*
      ** Wait for end of label.
      */
      if( comment[i] == '\n' ||
          comment[i] == '\0' )
      {
        comment[i] = '\0';
        ++n;
        istate = 0;
      }
    }
  }

  *numEntry = n;

  return( ce );
}


#ifdef ANSI_C
char *findSignal( char *src, char target )
#else
char *findSignal( src, target )
char *src;
char target;
#endif
{
  int i, j;
  int len;
  static char substring[1024];

  substring[0] = '\0';
  len = strlen( src );

  /*
  ** Search for target character.
  */
  for( i = 0; i < len; ++i )
  {
    if( src[i] == target )
    {
      j = i + 1;
      break;
    }
    if( i == len - 1 )
    {
      return( NULL );
    }
  }

  /*
  ** Skip spaces and colon.
  */
  for( i = j; i < len; ++i )
  {
    if( src[i] != ':' &&
        src[i] != '=' &&
        src[i] != ' ' )
    {
      j = i;
      break;
    }
    if( i == len - 1 )
    {
      return( NULL );
    }
  }

  /*
  ** Copy value.
  */
  for( i = j; i < len; ++i )
  {
    if( src[i] == ' ' ||
        src[i] == ',' ||
        src[i] == '\0' ||
        i - j == 1023 )
    {
      substring[i-j] = '\0';
      break;
    }
    substring[i-j] = src[i];
  }

  return( substring );
}



