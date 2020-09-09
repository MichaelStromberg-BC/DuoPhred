/** writeSCF3.c **/

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
**    * writeSCF3.c                                                *         **
**    * benefits from ideas in code written by LaDeana Hillier and *         **
**    * Tim Gleeson.                                               *         **
**                                                                           **
*******************************************************************************
*/






#define VERSION_NO "3.00"

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
  uint4 private_size;          /* private data block */
  uint4 private_offset;
  uint4 spare[18];             /* Unused */
} Header;

#define CSET_DEFAULT 0         /* {A,C,G,T,-} */
#define CSET_STADEN  1
#define CSET_NC_IUB  2
#define CSET_ALF     3         /* extended NC_IUB */
#define CSET_ABI     4         /* {A,C,G,T,N} */

/*
** Type definition for the comments
*/
typedef char Comments;

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
** writeSCF3
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
int writeSCF3( char *fn, int prec, PhredData *phredData )
#else
int writeSCF3( fn, prec, phredData )
char *fn;
int prec;
PhredData *phredData;
#endif
{

  int i, j;
  int scaleFlag;
  FLOAT min, max;
  FLOAT scale_offset;
  FLOAT scale_factor;
  int numByte;
  int numBase;
  int numPoint;
  int begWrite;
  int numWrite;
  uint4 *ui4;
  uint1 *ui1;
  char cmntString[PHRED_MAX_STRING_LEN];
  FILE *fp;
  Header header;
  Comments comments[12*PHRED_MAX_STRING_LEN];
  Option *option;

  option = getOption();

  /*
  ** Open output file.
  */
  if( ( fp = fopen( fn, "wb+" ) ) == NULL ) 
  {
    fprintf( stderr, "writeSCF: unable to open file %s\n", fn );
    return( ERROR );
  }

  /*
  ** Copy some often used values.
  */
  numBase  = phredData->numBase[phredData->dataSet[0]];
  numPoint = phredData->numPoint;

  pstrcpy( cmntString, phredData->comment, PHRED_MAX_STRING_LEN );

  if( option->trimSCFData == 1 )
  {
    pstrcat( cmntString, " trimmed", PHRED_MAX_STRING_LEN );
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
  sprintf( comments, "SIGN=A=%d,C=%d,G=%d,T=%d\nSPAC=%6.2f\nPRIM=%d\nMACH=%.*s\nDYEP=%.*s\nNAME=%.*s\nLANE=%d\nGELN=%.*s\nPROC=%.*s\nRTRK=%.*s\nCONV=%.16s%.*s\nCOMM=%.*s\nSRCE=%.*s",
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

  /*
  ** We may need to trim `base' information.
  */
  if( option->trimSCFData == 1 )
  {
    begWrite = phredData->leftTrimPoint;
    numWrite = numBase - phredData->rghtTrimPoint - phredData->leftTrimPoint;
  }
  else
  {
    begWrite = 0;
    numWrite= numBase;
  }

  /*
  ** Initialize header structure.
  */
  header.magic_number     = SCF_MAGIC;
  header.samples          = numPoint;
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
  header.bases_offset = (uint4)( header.samples_offset +
                                 4 * header.samples *
                                 ( ( prec == 2 ) ?
                                   sizeof( uint2 ) :
                                   sizeof( uint1 ) ) );

  header.comments_size    = (uint4)strlen( comments ) + 1;
  header.comments_offset  = (uint4)( header.bases_offset +
                                     header.bases * 12 );
  strncpy( header.version, VERSION_NO , 4 );
  header.sample_size      = prec;
  header.code_set         = CSET_DEFAULT;
  header.private_size     = 0;
  header.private_offset   = (uint4)( header.comments_offset +
                                     header.comments_size );

  for( i = 0; i < 18; ++i )
  {
    header.spare[i] = 0;
  }

  /*
  ** Swap bytes.
  */
  header.magic_number     = outSwpUint4( header.magic_number );
  header.samples          = outSwpUint4( header.samples );
  header.samples_offset   = outSwpUint4( header.samples_offset );
  header.bases            = outSwpUint4( header.bases );
  header.bases_left_clip  = outSwpUint4( header.bases_left_clip );
  header.bases_right_clip = outSwpUint4( header.bases_right_clip );
  header.bases_offset     = outSwpUint4( header.bases_offset );
  header.comments_size    = outSwpUint4( header.comments_size );
  header.comments_offset  = outSwpUint4( header.comments_offset );
  header.sample_size      = outSwpUint4( header.sample_size );
  header.code_set         = outSwpUint4( header.code_set );
  header.private_size     = outSwpUint4( header.private_size );
  header.private_offset   = outSwpUint4( header.private_offset );
  for( i = 0; i < 18; ++i )
  {
    header.spare[i] = outSwpUint4( header.spare[i] );
  }

  /*
  ** Write header.
  */
  if( fwrite( &header, sizeof( Header ), 1, fp ) != 1 )
  {
    fprintf( stderr,
             "writeSCF3: error: unable to write header\n" );
    fclose( fp );
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
  if( prec != 1 && phredData->maxTraceValue > 65535.0 )
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
      scale_factor = 65535.0 / ( max - min );
    }
  }

  if( option->verboseOption == 1 && option->verboseLevel > 60 )
  {
    fprintf( stderr,
             "writeSCF3: SCF version:     3.0\n" );
    fprintf( stderr,
             "writeSCF3: trace precision: %d\n",
             prec );
    fprintf( stderr,
             "writeSCF3: scaling:         %s\n",
             scaleFlag == 0 ? "no" : "yes" );
    if( scaleFlag == 1 )
    {
      fprintf( stderr,
               "writeSCF3: scale offset:    %f\n",
               scale_offset );
      fprintf( stderr,
               "writeSCF3: scale factor:    %f\n",
               scale_factor );
    }
    if( option->trimSCFData == 1 )
    {
      fprintf( stderr,
               "writeSCF3: trim bases: retain %d to %d (first base is zero)\n",
               phredData->leftTrimPoint,
               numBase - phredData->rghtTrimPoint - 1 );
    }
  }

#define SCALE( V, F, O )	( ( (V) - (O) ) * (F) )

  /*
  ** Write trace data.
  */
  if( numPoint > 0 )
  {
    if( prec == 1 )
    {
      uint1 prev, itmp;
      uint1 *itrace;

      /*
      ** Allocate temporary storage for trace.
      */
      numByte = phredData->numPoint * sizeof( uint1 );
      itrace = (uint1 *)malloc( numByte );
      if( itrace == NULL )
      {
        fprintf( stderr, "writeSCF3: error: unable to allocate memory\n" );
        fclose( fp );
        return( ERROR );
      }

      if( scaleFlag == 1 )
      {
        for( j = 0; j < 4; ++j )
        {
          prev = 0;
          for( i = 0; i < phredData->numPoint; ++i )
          {
            itmp = (uint1)SCALE( phredData->trace[j][i],
                                 scale_factor,
                                 scale_offset );
            itrace[i] = itmp - prev;
            prev = itmp;
          }
          prev = 0;
          for( i = 0; i < phredData->numPoint; ++i )
          {
            itmp = itrace[i];
            itrace[i] = itmp - prev;
            prev = itmp;
          }
          if( fwrite( itrace, sizeof( uint1 ), phredData->numPoint, fp ) != phredData->numPoint )
          {
            fprintf( stderr, "writeSCF3: error: unable to write trace\n" );
            free( itrace );
            fclose( fp );
            return( ERROR );
          }
        }
      }
      else
      {
        for( j = 0; j < 4; ++j )
        {
          prev = 0;
          for( i = 0; i < phredData->numPoint; ++i )
          {
            itmp = (uint1)phredData->trace[j][i];
            itrace[i] = itmp - prev;
            prev = itmp;
          }
          prev = 0;
          for( i = 0; i < phredData->numPoint; ++i )
          {
            itmp = itrace[i];
            itrace[i] = itmp - prev;
            prev = itmp;
          }
          if( fwrite( itrace, sizeof( uint1 ), phredData->numPoint, fp ) != phredData->numPoint )
          {
            fprintf( stderr, "writeSCF3: error: unable to write trace\n" );
            free( itrace );
            fclose( fp );
            return( ERROR );
          }
        }
      }

      free( itrace );
    }
    else
    {
      uint2 prev, itmp;
      uint2 *itrace;

      /*
      ** Allocate temporary storage for trace.
      */
      numByte = phredData->numPoint * sizeof( uint2 );
      itrace = (uint2 *)malloc( numByte );
      if( itrace == NULL )
      {
        fprintf( stderr, "writeSCF3: error: unable to allocate memory\n" );
        fclose( fp );
        return( ERROR );
      }

      if( scaleFlag == 1 )
      {
        for( j = 0; j < 4; ++j )
        {
          prev = 0;
          for( i = 0; i < phredData->numPoint; ++i )
          {
            itmp = (uint2)SCALE( phredData->trace[j][i],
                                 scale_factor,
                                 scale_offset );
            itrace[i] = itmp - prev;
            prev = itmp;
          }
          prev = 0;
          for( i = 0; i < phredData->numPoint; ++i )
          {
            itmp = itrace[i];
            itrace[i] = outSwpUint2( itmp - prev );
            prev = itmp;
          }
          if( fwrite( itrace, sizeof( uint2 ), phredData->numPoint, fp ) != phredData->numPoint )
          {
            fprintf( stderr, "writeSCF3: error: unable to write trace\n" );
            free( itrace );
            fclose( fp );
            return( ERROR );
          }
        }
      }
      else
      {
        for( j = 0; j < 4; ++j )
        {
          prev = 0;
          for( i = 0; i < phredData->numPoint; ++i )
          {
            itmp = (uint2)phredData->trace[j][i];
            itrace[i] = itmp - prev;
            prev = itmp;
          }
          prev = 0;
          for( i = 0; i < phredData->numPoint; ++i )
          {
            itmp = itrace[i];
            itrace[i] = outSwpUint2( itmp - prev );
            prev = itmp;
          }
          if( fwrite( itrace, sizeof( uint2 ), phredData->numPoint, fp ) != phredData->numPoint )
          {
            fprintf( stderr, "writeSCF3: error: unable to write trace\n" );
            free( itrace );
            fclose( fp );
            return( ERROR );
          }
        }
      }

      free( itrace );
    }
  }

  /*
  ** Write base locations.
  */
  if( numBase > 0 )
  {
    numByte = numBase * sizeof( uint4 );
    ui4 = (uint4 *)malloc( numByte );
    if( ui4 == NULL )
    {
      fprintf( stderr, "writeSCF3: error: unable to allocate memory\n" );
      fclose( fp );
    }

    for( i = 0; i < numBase; ++i )
    {
      ui4[i] = outSwpUint4( (uint4)phredData->baseLoc[phredData->dataSet[0]][i] );
    }
    if( fwrite( &(ui4[begWrite]), sizeof( uint4 ), numWrite, fp ) != numWrite )
    {
      fprintf( stderr, "writeSCF3: error: unable to write base locations\n" );
      free( ui4 );
      fclose( fp );
      return( ERROR );
    }
    free( ui4 );
  }

  /*
  ** Write quality values.
  */
  if( numBase > 0 )
  {
    /*
    ** Allocate memory for storing quality values.
    */
    numByte = numBase * sizeof( uint1 );
    ui1 = (uint1 *)malloc( numByte );
    if( ui1 == NULL )
    {
      fprintf( stderr, "writeSCF3: error: unable to allocate memory\n" );
      fclose( fp );
    }

    /*
    ** Store quality values.
    */
    numByte = numBase * sizeof( uint1 );
    memset( ui1, 0, numByte );
    for( i = 0; i < numBase; ++i )
    {
      if( phredData->base[phredData->dataSet[0]][i] == 'A' ||
          phredData->base[phredData->dataSet[0]][i] == 'N' )
      {
        ui1[i] = phredData->baseQual[phredData->dataSet[0]][i];
      }
    }
    if( fwrite( &(ui1[begWrite]), sizeof( uint1 ), numWrite, fp ) != numWrite )
    {
      fprintf( stderr, "writeSCF3: error: unable to write A base probabilities\n" );
      free( ui1 );
      fclose( fp );
      return( ERROR );
    }

    numByte = numBase * sizeof( uint1 );
    memset( ui1, 0, numByte );
    for( i = 0; i < numBase; ++i )
    {
      if( phredData->base[phredData->dataSet[0]][i] == 'C' ||
          phredData->base[phredData->dataSet[0]][i] == 'N' )
      {
        ui1[i] = phredData->baseQual[phredData->dataSet[0]][i];
      }
    }
    if( fwrite( &(ui1[begWrite]), sizeof( uint1 ), numWrite, fp ) != numWrite )
    {
      fprintf( stderr, "writeSCF3: error: unable to write C base probabilities\n" );
      free( ui1 );
      fclose( fp );
      return( ERROR );
    }

    numByte = numBase * sizeof( uint1 );
    memset( ui1, 0, numByte );
    for( i = 0; i < numBase; ++i )
    {
      if( phredData->base[phredData->dataSet[0]][i] == 'G' ||
          phredData->base[phredData->dataSet[0]][i] == 'N' )
      {
        ui1[i] = phredData->baseQual[phredData->dataSet[0]][i];
      }
    }
    if( fwrite( &(ui1[begWrite]), sizeof( uint1 ), numWrite, fp ) != numWrite )
    {
      fprintf( stderr, "writeSCF3: error: unable to write G base probabilities\n" );
      free( ui1 );
      fclose( fp );
      return( ERROR );
    }

    numByte = numBase * sizeof( uint1 );
    memset( ui1, 0, numByte );
    for( i = 0; i < numBase; ++i )
    {
      if( phredData->base[phredData->dataSet[0]][i] == 'T' ||
          phredData->base[phredData->dataSet[0]][i] == 'N' )
      {
        ui1[i] = phredData->baseQual[phredData->dataSet[0]][i];
      }
    }
    if( fwrite( &(ui1[begWrite]), sizeof( uint1 ), numWrite, fp ) != numWrite )
    {
      fprintf( stderr, "writeSCF3: error: unable to write T base probabilities\n" );
      free( ui1 );
      fclose( fp );
      return( ERROR );
    }
  }

  /*
  ** Write bases.
  */
  if( numBase > 0 )
  {
    if( fwrite( &(phredData->base[phredData->dataSet[0]][begWrite]),
                sizeof( char ),
                numWrite,
                fp ) != numWrite )
    {
      fprintf( stderr, "writeSCF3: error: unable to write bases\n" );
      free( ui1 );
      fclose( fp );
      return( ERROR );
    }
  }

  /*
  ** Write spare.
  */
  if( numBase > 0 )
  {
    numByte = numBase * sizeof( uint1 );
    memset( ui1, 0, numByte );
    if( fwrite( &(ui1[begWrite]), sizeof( uint1 ), numWrite, fp ) != numWrite )
    {
      fprintf( stderr, "writeSCF3: error: unable to write base spare values\n" );
      free( ui1 );
      fclose( fp );
      return( ERROR );
    }
    if( fwrite( &(ui1[begWrite]), sizeof( uint1 ), numWrite, fp ) != numWrite )
    {
      fprintf( stderr, "writeSCF3: error: unable to write base spare values\n" );
      free( ui1 );
      fclose( fp );
      return( ERROR );
    }
    if( fwrite( &(ui1[begWrite]), sizeof( uint1 ), numWrite, fp ) != numWrite )
    {
      fprintf( stderr, "writeSCF3: error: unable to write base spare values\n" );
      free( ui1 );
      fclose( fp );
      return( ERROR );
    }

    /*
    ** Free memory.
    */
    free( ui1 );
  }

  /*
  ** Write comments.
  */
  if( write_scf_comment( fp, comments, (size_t)( strlen( comments ) + 1 ) ) == ERROR )
  {
    fprintf( stderr, "writeSCF3: error: unable to write comments\n" );
    return( ERROR );
  }

  /*
  ** Done...prepare to return.
  */
  fclose( fp );
    
  return( OK );
}


