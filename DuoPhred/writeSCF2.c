/** writeSCF2.c **/

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
**    * writeSCF2.c                                                *         **
**    * benefits from ideas in code written by LaDeana Hillier and *         **
**    * Tim Gleeson.                                               *         **
**                                                                           **
*******************************************************************************
*/

/*
** Note: try to have phred write a 'empty' output files
**       on error.
*/

#define VERSION_NO "2.00"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "phred.h"
#include "rwUtil.h"

#define SCF_MAGIC (((((uint4)'.'<<8)+(uint4)'s'<<8)+(uint4)'c'<<8)+(uint4)'f')

/*
** Type definition for the Header structure
*/
typedef struct
{
  uint4 magic_number;          /* SCF_MAGIC */
  uint4 samples;               /* Number of elements in Samples matrix */
  uint4 samples_offset;        /* Byte offset from start of file */
  uint4 bases;                 /* Number of bases in Bases matrix */
  uint4 bases_left_clip;       /* Number of bases in left clip (vector)*/
  uint4 bases_right_clip;      /* Number of bases in right clip (unreliable) */
  uint4 bases_offset;          /* Byte offset from start of file */
  uint4 comments_size;         /* Number of bytes in Comment section */
  uint4 comments_offset;       /* Byte offset from start of file */
  char version[4];             /* "version.revision" */
  uint4 sample_size;           /* precision of samples (in bytes) */
  uint4 code_set;              /* uncertainty codes used */
  uint4 spare[20];             /* Unused */
} Header;

#define CSET_DEFAULT 0         /* {A,C,G,T,-} */
#define CSET_STADEN  1
#define CSET_NC_IUB  2
#define CSET_ALF     3         /* extended NC_IUB */
#define CSET_ABI     4         /* {A,C,G,T,N} */

/*
** Type definition for the Sample data
*/
typedef unsigned char byte;

typedef struct
{
  byte sample_A;               /* Sample for A trace */
  byte sample_C;               /* Sample for C trace */
  byte sample_G;               /* Sample for G trace */
  byte sample_T;               /* Sample for T trace */
} Samples1;

typedef struct
{
  uint2 sample_A;     /* Sample for A trace */
  uint2 sample_C;     /* Sample for C trace */
  uint2 sample_G;     /* Sample for G trace */
  uint2 sample_T;     /* Sample for T trace */
} Samples2;

/*
** Type definition for the sequence data
*/
typedef struct
{
  uint4 peak_index;            /* Index into Samples matrix for base position */
  uint1 prob_A;                /* Probability of it being an A */
  uint1 prob_C;                /* Probability of it being an C */
  uint1 prob_G;                /* Probability of it being an G */
  uint1 prob_T;                /* Probability of it being an T */
  char base;                   /* Base called */
  uint1 spare[3];              /* Spare */
} Bases;


/*
** Type definition for the comments
*/
typedef char Comments;

#ifdef ANSI_C
static int write_scf_header( FILE *fp, Header *h )
#else
static int write_scf_header( fp, h )
FILE *fp;
Header *h;
#endif
{
  int i;
  if( writeUint4( fp, h->magic_number )      == ERROR ) return( ERROR );
  if( writeUint4( fp, h->samples )           == ERROR ) return( ERROR );
  if( writeUint4( fp, h->samples_offset )    == ERROR ) return( ERROR );
  if( writeUint4( fp, h->bases )             == ERROR ) return( ERROR );
  if( writeUint4( fp, h->bases_left_clip )   == ERROR ) return( ERROR );
  if( writeUint4( fp, h->bases_right_clip )  == ERROR ) return( ERROR );
  if( writeUint4( fp, h->bases_offset )      == ERROR ) return( ERROR );
  if( writeUint4( fp, h->comments_size )     == ERROR ) return( ERROR );
  if( writeUint4( fp, h->comments_offset )   == ERROR ) return( ERROR );
  if( fwrite( h->version, sizeof( h->version ), 1, fp ) !=  1 )
  {
    return( ERROR );
  }
  if( writeUint4( fp, h->sample_size )       == ERROR ) return( ERROR );
  if( writeUint4( fp, h->code_set )          == ERROR ) return( ERROR );
  for( i = 0; i < 20; ++i )
  {
    if( writeUint4( fp, h->spare[i] ) == ERROR )
    {
      return( ERROR );
    }
  }
    return( OK );
}


#ifdef ANSI_C
static int write_scf_sample1( FILE *fp, Samples1 *s )
#else
static int write_scf_sample1( fp, s )
FILE *fp;
Samples1 *s;
#endif
{
  if( writeUint1( fp, s->sample_A ) == ERROR ) return( ERROR );
  if( writeUint1( fp, s->sample_C ) == ERROR ) return( ERROR );
  if( writeUint1( fp, s->sample_G ) == ERROR ) return( ERROR );
  if( writeUint1( fp, s->sample_T ) == ERROR ) return( ERROR );

  return( OK );
}


#ifdef ANSI_C
static int write_scf_sample2( FILE *fp, Samples2 *s )
#else
static int write_scf_sample2( fp, s )
FILE *fp;
Samples2 *s;
#endif
{
  if( writeUint2( fp, s->sample_A ) == ERROR ) return( ERROR );
  if( writeUint2( fp, s->sample_C ) == ERROR ) return( ERROR );
  if( writeUint2( fp, s->sample_G ) == ERROR ) return( ERROR );
  if( writeUint2( fp, s->sample_T ) == ERROR ) return( ERROR );

  return( OK );
}





#ifdef ANSI_C
static int write_scf_base( FILE *fp, Bases *b )
#else
static int write_scf_base( fp, b )
FILE *fp;
Bases *b;
#endif
{
  if( writeUint4( fp, b->peak_index )    == ERROR ) return( ERROR );
  if( writeUint1( fp, b->prob_A )        == ERROR ) return( ERROR );
  if( writeUint1( fp, b->prob_C )        == ERROR ) return( ERROR );
  if( writeUint1( fp, b->prob_G )        == ERROR ) return( ERROR );
  if( writeUint1( fp, b->prob_T )        == ERROR ) return( ERROR );
  if( writeUint1( fp, (uint1)b->base)    == ERROR ) return( ERROR );
  if( writeUint1( fp, b->spare[0] )      == ERROR ) return( ERROR );
  if( writeUint1( fp, b->spare[1] )      == ERROR ) return( ERROR );
  if( writeUint1( fp, b->spare[2] )      == ERROR ) return( ERROR );

  return( OK );
}

#ifdef ANSI_C
static int write_scf_comment( FILE *fp, Comments *c, size_t l )
#else
static int write_scf_comment( fp, c, l )
FILE *fp;
Comments *c;
size_t l;
#endif
{
  if( fwrite( c, l, 1, fp ) != 1 )
  {
    return( ERROR );
  }
  return( OK );
}

/*
** writeSCF2
**
** Purpose
** =======
** write SCF file from Seq structure
**
** Arguments
** =========
**
**         fn:             (input)  char *
**                         Pointer to string of output
**                         file name.
**
**         prec            (input)  int
**                         Trace data sample size.
**                         1=1 byte  2=2 byte
**
**         seq             (input)  seq
**                         Pointer to sequence structure.
**
*/
#ifdef ANSI_C
int writeSCF2( char *fn, int prec, PhredData *phredData )
#else
int writeSCF2( fn, prec, phredData )
char *fn;
int prec;
PhredData *phredData;
#endif
{

  int i;
  int numBase;
  int scaleFlag;
  int begWrite;
  int endWrite;
  int numWrite;
  FLOAT min, max;
  FLOAT scale_offset;
  FLOAT scale_factor;
  char cmntString[PHRED_MAX_STRING_LEN];
  FILE *fp;
  Header header;
  Bases base;
  Comments comments[12*PHRED_MAX_STRING_LEN];
  Option *option;

  option = getOption();

  /*
  ** Open output file.
  */
  if( ( fp = fopen( fn, "wb+" ) ) == NULL ) 
  {
    fprintf( stderr, "writeSCF2: unable to open file %s\n", fn );
    return( ERROR );
  }

  /*
  ** Prepare comment string for flagging `trimmed'
  ** sequence.
  */
  pstrcpy( cmntString, phredData->comment, PHRED_MAX_STRING_LEN );

  if( cmntString[strlen( cmntString )-1] == '\n' )
  {
    cmntString[strlen( cmntString )-1] = '\0';
  }

  if( option->trimSCFData == 1 )
  {
    if( strlen( cmntString ) > 0 )
    {
      pstrcat( cmntString, " trimmed",  PHRED_MAX_STRING_LEN );
    }
    else
    {
      pstrcat( cmntString, "trimmed",  PHRED_MAX_STRING_LEN );
    }
  }

  /*
  ** Construct comment line(s).
  */
  sprintf( comments,
           "SIGN=A=%d,C=%d,G=%d,T=%d\nSPAC=%6.2f\nPRIM=%d\nMACH=%.*s\nDYEP=%.*s\nNAME=%.*s\nLANE=%2d\nGELN=%.*s\nPROC=%.*s\nRTRK=%.*s\nCONV=%.16s%.*s\nCOMM=%.*s\nSRCE=%.*s",
           phredData->signalStrength[0],
           phredData->signalStrength[1],
           phredData->signalStrength[2],
           phredData->signalStrength[3],
           phredData->avgSpacing,
           phredData->primerLoc,
           PHRED_MAX_STRING_LEN,
           phredData->machineName,
           PHRED_MAX_STRING_LEN,
           phredData->primerID,
           PHRED_MAX_STRING_LEN,
           phredData->sampleName,
           phredData->laneNumber,
           PHRED_MAX_STRING_LEN,
           phredData->gelName,
           PHRED_MAX_STRING_LEN,
           phredData->processing,
           PHRED_MAX_STRING_LEN,
           phredData->reTracker,
           "phred version=",
           PHRED_MAX_STRING_LEN,
           getVersion(),
           PHRED_MAX_STRING_LEN,
           cmntString,
           PHRED_MAX_STRING_LEN,
           phredData->source );

  numBase = phredData->numBase[phredData->dataSet[0]];

  /*
  ** Trim `bases' if necessary.
  */
  if( option->trimSCFData == 1 )
  {
    begWrite = phredData->leftTrimPoint;
    endWrite = numBase - phredData->rghtTrimPoint;
    numWrite = numBase - phredData->rghtTrimPoint - phredData->leftTrimPoint;
  }
  else
  {
    begWrite = 0;
    endWrite = numBase;
    numWrite = numBase;
  }
	
  /*
  ** Initialize header structure.
  */
  header.magic_number     = SCF_MAGIC;
  header.samples          = phredData->numPoint;
  header.samples_offset   = (uint4)sizeof( Header );

  if( option->trimSCFData == 1 )
  {
    header.bases            = numWrite;
    header.bases_left_clip  = 0;
    header.bases_right_clip = 0;
  }
  else
  {
    header.bases            = numBase;
    header.bases_left_clip  = phredData->leftTrimPoint;
    header.bases_right_clip = phredData->rghtTrimPoint;
  }
  header.bases_offset     = (uint4)( header.samples_offset +
                                     header.samples *
                                     ( ( prec == 2 ) ?
                                       sizeof( Samples2 ) :
                                       sizeof( Samples1 ) ) );

  header.comments_size    = (uint4)strlen( comments ) + 1;
  header.comments_offset  = (uint4)( header.bases_offset +
                                     header.bases *
                                     sizeof( Bases ) );
  strncpy( header.version, VERSION_NO , 4 );
  header.sample_size      = prec;
  header.code_set         = CSET_DEFAULT;

  for( i = 0; i < 20; ++i )
  {
    header.spare[i] = 0;
  }

  /*
  ** Write header.
  */
  if( write_scf_header( fp, &header ) == ERROR )
  {
    fclose( fp );
    fprintf( stderr,
             "writeSCF2: error: unable to write header\n" );
    return( ERROR );
  }

  /*
  ** Decide whether or not to scale traces.
  */
  scaleFlag = 0;
  if( option->scaleSCFTraceOption == 1 )
  {
    scaleFlag = 1;
  }
  else
  if( phredData->fileType == MD1Format ||
      phredData->fileType == MD2Format )
  {
    scaleFlag = 1;
  }
  else
  if( prec == 1 && phredData->maxTraceValue > 255.0 )
  {
    scaleFlag = 1;
  }
  else
  if( prec == 2 && phredData->maxTraceValue > 65535.0 )
  {
    scaleFlag = 1;
  }
  else
  if( prec != 1 && prec != 2 )
  {
    scaleFlag = 1;
  }


  /*
  ** Determine trace scaling values.
  */
  if( scaleFlag == 1 )
  {
    min = phredData->minTraceValue;
    max = phredData->maxTraceValue;

    scale_offset = min;

    if( min == max )
    {
      /*
      ** Avoid a divide by zero condition.
      ** The scale_offset sets trace to all zeroes
      ** in this case.
      */
      scale_factor = 1.0;
    }
    else
    if( prec == 1 )
    {
      scale_factor = 255.0 / ( max - min );
    }
    else
    if( prec == 2 )
    {
      scale_factor = 65535.0 / ( max - min );
    }
    else
    {
      scale_factor = 255.0 / ( max - min );
    }
  }

  if( option->verboseOption == 1 && option->verboseLevel > 60 )
  {
    fprintf( stderr,
             "writeSCF2: SCF version:     2.0\n" );
    fprintf( stderr,
             "writeSCF2: trace precision: %d\n",
             prec );
    fprintf( stderr,
             "writeSCF2: scaling:         %s\n",
             scaleFlag == 0 ? "no" : "yes" );
    if( scaleFlag == 1 )
    {
      fprintf( stderr,
               "writeSCF2: scale offset:    %f\n",
               scale_offset );
      fprintf( stderr,
               "writeSCF2: scale factor:    %f\n",
               scale_factor );
    }
    if( option->trimSCFData == 1 )
    {
      fprintf( stderr,
               "writeSCF2: trim bases: retain %d to %d (first base is zero)\n",
               phredData->leftTrimPoint,
               numBase - phredData->rghtTrimPoint - 1 );
    }
  }

#define SCALE( V, F, O )	( ( (V) - (O) ) * (F) )

  /*
  ** Write trace data.
  */
  switch( prec )
  {
    case 1:
    {
      Samples1 sample;
	
      for( i = 0; i < (int)header.samples; ++i )
      {
        if( scaleFlag == 1 )
        {
          sample.sample_A = (uint1)SCALE( phredData->trace[0][i],
                                          scale_factor,
                                          scale_offset );
          sample.sample_C = (uint1)SCALE( phredData->trace[1][i],
                                          scale_factor,
                                          scale_offset );
          sample.sample_G = (uint1)SCALE( phredData->trace[2][i],
                                          scale_factor,
                                          scale_offset );
          sample.sample_T = (uint1)SCALE( phredData->trace[3][i],
                                          scale_factor,
                                          scale_offset );
        }
        else
        {
          sample.sample_A = (uint1)phredData->trace[0][i];
          sample.sample_C = (uint1)phredData->trace[1][i];
          sample.sample_G = (uint1)phredData->trace[2][i];
          sample.sample_T = (uint1)phredData->trace[3][i];
        }
        if( write_scf_sample1( fp, &sample ) == ERROR )
        {
          fprintf( stderr,
                   "writeSCF2: error: unable to write trace\n" );
          fclose( fp );
          return( ERROR );
        }
      }
    }

    break;

    case 2:
    {
      Samples2 sample;

      for( i = 0; i < (int)header.samples; ++i )
      {
        if( scaleFlag == 1 )
        {
          sample.sample_A = (uint2)SCALE( phredData->trace[0][i],
                                          scale_factor,
                                          scale_offset );
          sample.sample_C = (uint2)SCALE( phredData->trace[1][i],
                                          scale_factor,
                                          scale_offset );
          sample.sample_G = (uint2)SCALE( phredData->trace[2][i],
                                          scale_factor,
                                          scale_offset );
          sample.sample_T = (uint2)SCALE( phredData->trace[3][i],
                                          scale_factor,
                                          scale_offset );
        }
        else
        {
          sample.sample_A = (uint2)phredData->trace[0][i];
          sample.sample_C = (uint2)phredData->trace[1][i];
          sample.sample_G = (uint2)phredData->trace[2][i];
          sample.sample_T = (uint2)phredData->trace[3][i];
        }
        if( write_scf_sample2( fp, &sample ) == ERROR )
        {
          fprintf( stderr,
                   "writeSCF2: error: unable to write trace\n" );
          fclose( fp );
          return( ERROR );
        }
      }
    }

    break;

    default:
    {
      Samples1 sample;

      for( i = 0; i < (int)header.samples; ++i )
      {
        if( scaleFlag == 1 )
        {
          sample.sample_A = (uint1)SCALE( phredData->trace[0][i],
                                          scale_factor,
                                          scale_offset );
          sample.sample_C = (uint1)SCALE( phredData->trace[1][i],
                                          scale_factor,
                                          scale_offset );
          sample.sample_G = (uint1)SCALE( phredData->trace[2][i],
                                          scale_factor,
                                          scale_offset );
          sample.sample_T = (uint1)SCALE( phredData->trace[3][i],
                                          scale_factor,
                                          scale_offset );
        }
        else
        {
          sample.sample_A = (uint1)phredData->trace[0][i];
          sample.sample_C = (uint1)phredData->trace[1][i];
          sample.sample_G = (uint1)phredData->trace[2][i];
          sample.sample_T = (uint1)phredData->trace[3][i];
        }
        if( write_scf_sample1( fp, &sample ) == ERROR )
        {
          fprintf( stderr,
                   "writeSCF2: error: unable to write trace\n" );
          fclose( fp );
          return( ERROR );
        }
      }
    }

    break;
  }

  /*
  ** Write bases.
  */
  for( i = begWrite; i < endWrite; ++i )
  {
    base.peak_index = phredData->baseLoc[phredData->dataSet[0]][i];
    base.base = phredData->base[phredData->dataSet[0]][i];
    base.spare[0] = 0;
    base.spare[1] = 0;
    base.spare[2] = 0;
    switch( base.base )
    {
      case 'A':
        base.prob_A = phredData->baseQual[phredData->dataSet[0]][i];
        base.prob_C = 0;
        base.prob_G = 0;
        base.prob_T = 0;
        break;

      case 'C':
        base.prob_A = 0;
        base.prob_C = phredData->baseQual[phredData->dataSet[0]][i];
        base.prob_G = 0;
        base.prob_T = 0;
        break;

      case 'G':
        base.prob_A = 0;
        base.prob_C = 0;
        base.prob_G = phredData->baseQual[phredData->dataSet[0]][i];
        base.prob_T = 0;
        break;

      case 'T':
        base.prob_A = 0;
        base.prob_C = 0;
        base.prob_G = 0;
        base.prob_T = phredData->baseQual[phredData->dataSet[0]][i];
        break;

      case 'N':
        base.prob_A = phredData->baseQual[phredData->dataSet[0]][i];
        base.prob_C = phredData->baseQual[phredData->dataSet[0]][i];
        base.prob_G = phredData->baseQual[phredData->dataSet[0]][i];
        base.prob_T = phredData->baseQual[phredData->dataSet[0]][i];
        break;

      default:
        break;
    }

    if( write_scf_base( fp, &base ) == ERROR )
    {
      fprintf( stderr,
               "writeSCF2: error: unable to write bases\n" );
      fclose( fp );
      return( ERROR );
    }
  }

  /*
  ** Write comments.
  */
  if( write_scf_comment( fp, comments, (size_t)header.comments_size ) == ERROR )
  {
    fprintf( stderr,
             "writeSCF2: error: unable to write comments\n" );
    fclose( fp );
    return( ERROR );
  }

  /*
  ** Done...prepare to return.
  */
  fclose( fp );
    
  return( OK );
}


