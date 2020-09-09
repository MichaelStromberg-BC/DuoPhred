/** rwUtil.c **/

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
**    * rwUtil.c                                                       *     **
**    * benefits from ideas in code written by LaDeana Hillier and     *     **
**    * Tim Gleeson.                                                   *     **
**                                                                           **
*******************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "typeDef.h"
#include "rwUtil.h"

#ifdef ANSI_C
int2 inSwpSint2( char *ptr )
#else
int2 inSwpSint2( ptr )
char *ptr;
#endif
{
  int2 i;

  i = *ptr & 0xff;
  i = i << 8;
  ++ptr;

  i |= *ptr & 0xff;

  return( i );
}

#ifdef ANSI_C
uint2 inSwpUint2( char *ptr )
#else
uint2 inSwpUint2( ptr )
char *ptr;
#endif
{
  uint2 i;

  i = *ptr & 0xff;
  i = i << 8;
  ++ptr;

  i |= *ptr & 0xff;

  return( i );
}

#ifdef ANSI_C
int4 inSwpSint4( char *ptr )
#else
int4 inSwpSint4( ptr )
char *ptr;
#endif
{
  int4 i;

  i = *ptr & 0xff;
  i = i << 8;
  ++ptr;

  i |= *ptr & 0xff;
  i = i << 8;
  ++ptr;

  i |= *ptr & 0xff;
  i = i << 8;
  ++ptr;

  i |= *ptr & 0xff;

  return( i );
}

#ifdef ANSI_C
uint4 inSwpUint4( char *ptr )
#else
uint4 inSwpUint4( ptr )
char *ptr;
#endif
{
  uint4 i;

  i = *ptr & 0xff;
  i = i << 8;
  ++ptr;

  i |= *ptr & 0xff;
  i = i << 8;
  ++ptr;

  i |= *ptr & 0xff;
  i = i << 8;
  ++ptr;

  i |= *ptr & 0xff;

  return( i );
}

#ifdef ANSI_C
int2 outSwpSint2( int2 i )
#else
int2 outSwpSint2( i )
int2 i;
#endif
{
  int2  n;
  char *ptr;

  ptr    = (char *) &n;
  ptr[0] = ( i >> 8 ) & 0xff;
  ptr[1] = ( i      ) & 0xff;

  return( n );
}

#ifdef ANSI_C
uint2 outSwpUint2( uint2 i )
#else
uint2 outSwpUint2( i )
uint2 i;
#endif
{
  uint2  n;
  char  *ptr;

  ptr    = (char *) &n;
  ptr[0] = ( i >> 8 ) & 0xff;
  ptr[1] = ( i      ) & 0xff;

  return( n );
}

#ifdef ANSI_C
int4 outSwpSint4( int4 i )
#else
int4 outSwpSint4( i )
int4 i;
#endif
{
  int4  n;
  char *ptr;

  ptr    = (char *) &n;
  ptr[0] = ( i >> 24 ) & 0xff;
  ptr[1] = ( i >> 16 ) & 0xff;
  ptr[2] = ( i >>  8 ) & 0xff;
  ptr[3] = ( i       ) & 0xff;

  return( n );
}

#ifdef ANSI_C
uint4 outSwpUint4( uint4 i )
#else
uint4 outSwpUint4( i )
uint4 i;
#endif
{
  uint4  n;
  char *ptr;

  ptr    = (char *) &n;
  ptr[0] = ( i >> 24 ) & 0xff;
  ptr[1] = ( i >> 16 ) & 0xff;
  ptr[2] = ( i >>  8 ) & 0xff;
  ptr[3] = ( i       ) & 0xff;

  return( n );
}

#ifdef ANSI_C
float int2float( int4 i )
#else
float int2float( i )
int4 i;
#endif
{
  int4  f;
  int4  e;
  int4  s;
  float v;

  f = i & ( ( 1 << 23 ) - 1 );
  e = ( i >> 23 ) & ( ( 1 << 8 ) - 1 );
  s = i >> 31;

  v = (float)( ( s ? -1.0 : 1.0 ) *
               exp( log( (double)2.0 ) * (double)( e - 150 ) ) *
               (double)( ( 1 << 23 ) + f ) );
  
  return( v );
}

#ifdef ANSI_C
int readSint1( FILE *fp, int1 *i )
#else
int readSint1( fp, i )
FILE *fp;
int1 *i;
#endif
{
  if( fread( (char *)i, sizeof( int1 ), 1, fp ) != 1 )
  {
    return( ERROR );
  }
  return( OK );
}

#ifdef ANSI_C
int readUint1( FILE *fp, uint1 *i )
#else
int readUint1( fp, i )
FILE *fp;
uint1 *i;
#endif
{
  if( fread( (char *)i, sizeof( uint1 ), 1, fp ) != 1 )
  {
    return( ERROR );
  }
  return( OK );
}

#ifdef ANSI_C
int readSint2( FILE *fp, int2 *i )
#else
int readSint2( fp, i )
FILE *fp;
int2 *i;
#endif
{
  char buf[sizeof(int2)];

  if( fread( buf, sizeof( int2 ), 1, fp ) != 1 )
  {
    return( ERROR );
  }

  *i = inSwpSint2( buf );

  return( OK );
}

#ifdef ANSI_C
int readUint2( FILE *fp, uint2 *i )
#else
int readUint2( fp, i )
FILE *fp;
uint2 *i;
#endif
{
  char buf[sizeof(uint2)];

  if( fread( buf, sizeof( uint2 ), 1, fp ) != 1 )
  {
    return( ERROR );
  }

  *i = inSwpUint2( buf );

  return( OK );
}

#ifdef ANSI_C
int readSint4( FILE *fp, int4 *i )
#else
int readSint4( fp, i )
FILE *fp;
int4 *i;
#endif
{
  char buf[sizeof(int4)];

  if( fread( buf, sizeof( int4 ), 1, fp ) != 1 )
  {
    return( ERROR );
  }

  *i = inSwpSint4( buf );

  return( OK );
}

#ifdef ANSI_C
int readUint4( FILE *fp, uint4 *i )
#else
int readUint4( fp, i )
FILE *fp;
uint4 *i;
#endif
{
  char buf[sizeof(uint4)];

  if( fread( buf, sizeof( uint4 ), 1, fp ) != 1 )
  {
    return( ERROR );
  }

  *i = inSwpUint4( buf );

  return( OK );
}

#ifdef ANSI_C
int writeSint1( FILE *fp, int1 i )
#else
int writeSint1( fp, i )
FILE *fp;
int1 i;
#endif
{
  if( fwrite( (char *) &i, sizeof( int1 ), 1, fp ) != 1 )
  {
    return( ERROR );
  }
  return( OK );
}

#ifdef ANSI_C
int writeUint1( FILE *fp, uint1 i )
#else
int writeUint1( fp, i )
FILE *fp;
uint1 i;
#endif
{
  if( fwrite( (char *) &i, sizeof( uint1 ), 1, fp ) != 1 )
  {
    return( ERROR );
  }
  return( OK );
}

#ifdef ANSI_C
int writeSint2( FILE *fp, int2 i )
#else
int writeSint2( fp, i )
FILE *fp;
int2 i;
#endif
{
  int2 n;

  n = outSwpSint2( i );

  if( fwrite( (char *) &n, sizeof( int2 ), 1, fp ) != 1 )
  {
    return( ERROR );
  }
  return( OK );
}

#ifdef ANSI_C
int writeUint2( FILE *fp, uint2 i )
#else
int writeUint2( fp, i )
FILE *fp;
uint2 i;
#endif
{
  uint2 n;

  n = outSwpUint2( i );

  if( fwrite( (char *) &n, sizeof( uint2 ), 1, fp ) != 1 )
  {
    return( ERROR );
  }
  return( OK );
}

#ifdef ANSI_C
int writeSint4( FILE *fp, int4 i )
#else
int writeSint4( fp, i )
FILE *fp;
int4 i;
#endif
{
  int4 n;

  n = outSwpSint4( i );

  if( fwrite( (char *) &n, sizeof( int4 ), 1, fp ) != 1 )
  {
    return( ERROR );
  }
  return( OK );
}

#ifdef ANSI_C
int writeUint4( FILE *fp, uint4 i )
#else
int writeUint4( fp, i )
FILE *fp;
uint4 i;
#endif
{
  uint4 n;

  n = outSwpUint4( i );

  if( fwrite( (char *) &n, sizeof( uint4 ), 1, fp ) != 1 )
  {
    return( ERROR );
  }
  return( OK );
}


