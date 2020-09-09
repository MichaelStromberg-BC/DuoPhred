/** rwUtil.h **/

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

#ifdef ANSI_C
int2  inSwpSint2( char *ptr );
int4  inSwpSint4( char *ptr );
uint2 inSwpUint2( char *ptr );
uint4 inSwpUint4( char *ptr );
int2  outSwpSint2( int2 i );
int4  outSwpSint4( int4 i );
uint2 outSwpUint2( uint2 i );
uint4 outSwpUint4( uint4 i );
float int2float( int4 i );
int readSint1( FILE *fp, int1 *i );
int readSint2( FILE *fp, int2 *i );
int readSint4( FILE *fp, int4 *i );
int readUint1( FILE *fp, uint1 *i );
int readUint2( FILE *fp, uint2 *i );
int readUint4( FILE *fp, uint4 *i );
int writeSint1( FILE *fp, int1 i );
int writeSint2( FILE *fp, int2 i );
int writeSint4( FILE *fp, int4 i );
int writeUint1( FILE *fp, uint1 i );
int writeUint2( FILE *fp, uint2 i );
int writeUint4( FILE *fp, uint4 i );
#else
int2  inSwpSint2();
int4  inSwpSint4();
uint2 inSwpUint2();
uint4 inSwpUint4();
int2  outSwpSint2();
int4  outSwpSint4();
uint2 outSwpUint2();
uint4 outSwpUint4();
float int2float();
int readSint1();
int readSint2();
int readSint4();
int readUint1();
int readUint2();
int readUint4();
int writeSint1();
int writeSint2();
int writeSint4();
int writeUint1();
int writeUint2();
int writeUint4();
#endif

