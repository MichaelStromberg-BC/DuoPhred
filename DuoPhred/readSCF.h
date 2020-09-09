/** readSCF.h **/

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

#ifndef READSCF_DEFINED
#define READSCF_DEFINED

#include    "chromatData.h"

/* ---- Constants ---- */
#define SCF_MAGIC (((((int4)'.'<<8)+(int4)'s'<<8)+(int4)'c'<<8)+(int4)'f')
#define scale(V,OLDMAX,NEWMAX) (int2)( (FLOAT) V * (FLOAT) NEWMAX / (FLOAT) OLDMAX )

/*
** Type definition for the Header structure.
*/
typedef struct
{
  uint4 magic_number;
  uint4 samples;
  uint4 samples_offset;
  uint4 bases;
  uint4 bases_left_clip;
  uint4 bases_right_clip;
  uint4 bases_offset;
  uint4 comments_size;
  uint4 comments_offset;
  char  version[4];
  uint4 sample_size;
  uint4 code_set;
  uint4 private_size;
  uint4 private_offset;
  uint4 spare[18];
} SCFHeader;

/*
** Type definition for the Sample data
*/
typedef struct
{
  uint1 sample_A;
  uint1 sample_C;
  uint1 sample_G;
  uint1 sample_T;
} Samples1;

typedef struct
{
  uint2 sample_A;
  uint2 sample_C;
  uint2 sample_G;
  uint2 sample_T;
} Samples2;

/*
** Type definition for the sequence data
*/
typedef struct
{
  uint4 peak_index;
  uint1 prob_A;
  uint1 prob_C;
  uint1 prob_G;
  uint1 prob_T;
  uint1 base;
  uint1 spare[3];
} Bases;


typedef struct
{
  char *label;
  char *value;
  int  used;
} CommentEntry;


#ifdef ANSI_C
ChromatData *readSCF( char *fn, int *status );
uint2 swpBytU2( char *ptr );
uint4 swpBytU4( char *ptr );
char *findValue( int len, char *comment, char *target, int *ipos );
char *findSignal( char *src, char target );
CommentEntry *parseSCFComment( int len, char *comment, int *numEntry );
int readSCF2( char *fn, FILE *fp, SCFHeader *header, ChromatData *chromatData, int numPoint, int numBase, int *status );
int readSCF3( char *fn, FILE *fp, SCFHeader *header, ChromatData *chromatData, int numPoint, int numBase, int *status );
#else
ChromatData *readSCF();
uint2 inSwpUint2();
uint4 inSwpUint4();
char *findValue();
char *findSignal();
CommentEntry *parseSCFComment();
int readSCF2();
int readSCF3();
#endif

#endif
