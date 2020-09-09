/** readABI.c **/

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
**    * readABI.c                                                      *     **
**    * benefits from ideas in code written by Alan Blanchard,         *     **
**    * LaDeana Hillier, and Tim Gleeson.                              *     **
**                                                                           **
*******************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <sys/types.h>
#include "phred.h"
#include "rwUtil.h"
#include "readABI.h"

#ifndef SEEK_SET
#define SEEK_SET	0
#endif

/*
** Store ABI data block index information.
*/
typedef struct
{
  int     occur;
  int     sernum;
  uint4  offset;
  uint4  sizwrd;
  uint4  numwrd;
  uint4  numbyt;
  char    *label;
} DataIndex;

/*
** Number of data block indices we need.
** The names of the data blocks are
** listed in "labels" below.  NUMINDEX
** must equal the number of entries in
** "labels".
*/
#define NUMINDEX	17

/*
** Array indices for the data blocks that
** interest us.
**/
#define TRACE1	0
#define TRACE2	1
#define TRACE3	2
#define TRACE4	3
#define BASMAP	4
#define BASES	5
#define BASPOS	6
#define SIGSTR	7
#define AVGSPC	8
#define PRIPOS	9
#define MCHNAM	10
#define DYEPRI	11
#define SMPNAM  12
#define THMPRT	13
#define LANENM  14
#define GELNAM  15
#define COMMNT  16

/*
** Protoypes.
*/
#ifdef ANSI_C
char *ourMalloc( int numByte );
int ourFree( char *ptr );
#else
char *ourMalloc();
int ourFree();
#endif

/*
** The data block index labels and
** serial numbers (occurrences) for
** each block.
*/
char *labels[] = { "DATA", "DATA", "DATA", "DATA",
                   "FWO_", "PBAS", "PLOC", "S/N%",
                   "SPAC", "PPOS", "MCHN", "PDMF",
                   "SMPL", "THUM", "LANE", "GELN",
                   "CMNT" };

int  blockSerNum[NUMINDEX] = { 9, 10, 11, 12,
                               1,  1,  1,  1,
                               1,  1,  1,  1,
                               1,  1,  1,  1,
                               1 };
                   
/*
** readIndex
**
** Purpose:
**
**        Read blocks that contain information about
**        the size and location of the data blocks
**        within the file.
**
** Parameters:
**
**        fp           pointer to open ABI chromat file
**
**        offset       number of bytes in the file "header"
**                     (0 for files lacking the Mac header
**                      128 for files with the Mac header.)
**
**        maxNumByte   maximum number of bytes to read
**
**
*/
#ifdef ANSI_C
static DataIndex *readIndex( FILE *fp, int offset, int *maxNumByte )
#else
static DataIndex *readIndex( fp, offset, maxNumByte )
FILE *fp;
int offset;
int *maxNumByte;
#endif
{
  int i;
  int serNum;
  int datType;
  int numIndex;
  int corruptFlag;
  int iBlk;
  uint4 fstBlk;
  uint2 sizBlk;
  uint4 numBlk;
  static DataIndex dataIndex[NUMINDEX];
  char *buffer;
  char *cptr;
  Option *option;

  option = getOption();

  /*
  ** Initialize data index structure.
  */
  for( i = 0; i < NUMINDEX; ++i )
  {
    dataIndex[i].occur  = -1;
    dataIndex[i].offset = 0;
    dataIndex[i].numbyt = 0;
    dataIndex[i].numwrd = 0;
    dataIndex[i].label = labels[i];
    dataIndex[i].sernum = blockSerNum[i];
  }

  /*
  ** Read block size and allocate buffer memory.
  */
  if( fseek( fp, offset + 16, SEEK_SET ) != 0 ||
      readUint2( fp, &sizBlk ) == ERROR )
  {
    return( NULL );
  }

  /*
  ** Read number of header blocks.
  */
  if( fseek( fp, offset + 18, SEEK_SET ) != 0 ||
      readUint4( fp, &numBlk ) == ERROR )
  {
    return( NULL );
  }

  if( ( buffer = (char *)malloc( (size_t)sizBlk ) ) == NULL )
  {
    fprintf( stderr,
             "readABI: unable to allocate memory: sizBlk: %d\n",
             sizBlk );
    return( NULL );
  }

  /*
  ** Read first block index offset.
  */
  if( fseek( fp, offset + 26, SEEK_SET ) != 0 ||
      readUint4( fp, &fstBlk ) == ERROR )
  {
    free( buffer );
    return( NULL );
  }
  fstBlk += offset;

  /*
  ** Position file pointer at first data block index.
  */
  if( fseek( fp, (size_t)fstBlk, SEEK_SET ) != 0 )
  {
    free( buffer );
    return( NULL );
  }

  /*
  ** Step down through the sizBlk byte data block indices and
  ** match block label with a "dataIndex[].label".  Store
  ** matches.
  */
  iBlk = 0;
  numIndex = 0;
  *maxNumByte = 0;
  corruptFlag = 0;
  while( iBlk < numBlk &&
         numIndex < NUMINDEX &&
         fread( buffer, (size_t)sizBlk, (size_t)1, fp ) == (size_t)1 )
  {

    /*
    ** Check for corrupted id.
    */
    if( buffer[0] < 32 || buffer[0] > 126 ||
        buffer[1] < 32 || buffer[1] > 126 ||
        buffer[2] < 32 || buffer[2] > 126 ||
        buffer[3] < 32 || buffer[3] > 126 )
    {
      corruptFlag = 1;
    }

    /*
    ** Back to the work at hand.
    */
    cptr = buffer + 4;
    serNum = inSwpUint4( cptr );
    for( i = 0; i < NUMINDEX; ++i )
    {
      if( strncmp( dataIndex[i].label, buffer, 4 ) == 0 &&
          dataIndex[i].sernum == serNum )
      {
        dataIndex[i].occur = 1;
        cptr = buffer + 8;
        datType = (int)inSwpUint2( cptr );
        cptr = buffer + 10;
        dataIndex[i].sizwrd = (uint4)inSwpUint2( cptr );
        cptr = buffer + 12;
        dataIndex[i].numwrd = inSwpUint4( cptr );
        cptr = buffer + 16;
        dataIndex[i].numbyt = inSwpUint4( cptr );
        cptr = buffer + 20;

        /*
        ** If the data is 4 or less bytes, the data is stored in
        ** offset, otherwise the offset contains the pointer to
        ** the data.  The pointer is a file offset.
        */
        if( dataIndex[i].numbyt <= 4 )
        {
          dataIndex[i].offset = 0;
          memcpy( &(dataIndex[i].offset), cptr, 4 );
        }
        else
        {
          dataIndex[i].offset = offset + inSwpUint4( cptr );
        }
        if( (int)dataIndex[i].numbyt > *maxNumByte )
        {
          *maxNumByte = dataIndex[i].numbyt;
        }
        ++numIndex;
        break;
      }
    }
    ++iBlk;
  }

  if( corruptFlag == 1 )
  {
    fprintf( stderr,
             "readABI: warning: probably corrupt file: skip file\n" );
    return( NULL );
  }

  free( buffer );
  return( dataIndex );
}


/* ---- Exports ---- */

/*
** status:
**
**   0 = OK
**   1 = file reading error
**   2 = no trace (and no bases assumed)
**   3 = no bases (but there is trace)
**  -1 = fatal error
*/

/*
** Read the ABI format sequence with name `fn' into `seq'.
*/
#ifdef ANSI_C
ChromatData *readABI( char *fn, int *status )
#else
ChromatData *readABI( fn, status )
char *fn;
int *status;
#endif
{
  int4 i;
  int occur;
  int maxNumByte;
  int crFlg;
  uint4 numPoint;
  uint4 numBase;
  uint4 numwrd;
  int baseChar[128];
  size_t numbyt;
  size_t offset;
  DataIndex *dataIndex;
  int4 hoff;
  long int lhoff[2];
  char filID[8];
  char baseMap[4];
  char *buffer;
  char *cptr;
  FLOAT *fptr;
  ChromatData *chromatData;
  FILE *fp;
  Option *option;

  option = getOption();

  baseChar['A'] = 0;
  baseChar['C'] = 1;
  baseChar['G'] = 2;
  baseChar['T'] = 3;

  lhoff[0] = 0L;
  lhoff[1] = 128L;

  *status = 0;

  crFlg = 0;

  /*
  ** Open file for reading.
  */
  fp = fopen( fn, "rb" );
  if( fp == NULL )
  {
    fprintf( stderr, "readABI: unable to open file %s\n", fn );
    *status = 1;
    return( NULL );
  }

  /*
  ** Does the file have the 128 byte header?
  */
  hoff = -1;
  for( i = 0; i < 2; ++i )
  {
    if( fseek( fp, lhoff[i], 0 ) != 0 ||
        !fread( filID, 1, 4, fp ) )
    {
      continue;
    }
    if( strncmp( filID, "ABIF", 4 ) == 0 )
    {
      hoff = (int4)lhoff[i];
      break;
    }
  }
  if( hoff == -1 )
  {
    fprintf( stderr, "readABI: unable to ID file %s\n", fn );
    fclose( fp );
    *status = 1;
    return( NULL );
  }

  if( option->verboseOption ==1 && option->verboseLevel >= 32 )
  {
    if( hoff == 0 )
    {
      fprintf( stderr,
               "readABI: ABI file lacks Mac header: %s\n",
               fn );
    }
    else
    if( hoff == 128 )
    {
      fprintf( stderr,
               "readABI: ABI file has Mac header: %s\n",
               fn );
    }
  }


  /*
  ** Read data block indices.
  */
  dataIndex = readIndex( fp, hoff, &maxNumByte );
  if( dataIndex == NULL || maxNumByte == 0 )
  {
    fprintf( stderr, "unable to read header\n" );
    fclose( fp );

    numPoint = 0;
    numBase  = 0;

    chromatData = allocChromatData( numPoint, numBase );
    if( chromatData == NULL )
    {
      fprintf( stderr,
               "readABI: unable to allocate memory: numPoint: %d\n",
               numPoint );
      fclose( fp );
      *status = -1;
      return( NULL );
    }

    chromatData->fileType           = ABIFormat;
    chromatData->primerLoc          = 0;
    chromatData->avgSpacing         = 0.0;
    chromatData->machineName[0]     = '\0';
    chromatData->sampleName[0]      = '\0';
    chromatData->primerID[0]        = '\0';
    chromatData->signalStrength[0]  = 0;
    chromatData->signalStrength[1]  = 0;
    chromatData->signalStrength[2]  = 0;
    chromatData->signalStrength[3]  = 0;
    chromatData->gelName[0]         = '\0';
    chromatData->laneNumber         = 0;
    chromatData->processing[0]      = '\0';
    chromatData->reTracker[0]       = '\0';
    chromatData->comment[0]         = '\0';
    chromatData->convProg[0]        = '\0';
    chromatData->thumbPrint[0]      = 0;
    pstrcpy( chromatData->fileName, fn, PHRED_PATH_MAX );
    strcpy( chromatData->source, "ABI 373A or 377" );
    memset( chromatData->thumbPrint, 0, 10 * sizeof( char ) );

    *status = 1;
    return( chromatData );
  }

  /*
  ** Allocate buffer.
  */
  buffer = (char *)ourMalloc( maxNumByte );
  if( buffer == NULL )
  {
    fprintf( stderr,
             "readABI: unable to allocate memory: maxNumByte: %d\n",
             maxNumByte );
    *status = -1;
    fclose( fp );
    return( NULL );
  }

  /*
  ** Get the numbers of trace points and bases.
  */
  numPoint = dataIndex[TRACE1].numwrd;
  numBase = dataIndex[BASES].numwrd;

  /*
  ** Allocate chromatData.
  */
  chromatData = allocChromatData( numPoint, numBase );
  if( chromatData == NULL )
  {
    fprintf( stderr, "readABI: unable to allocate memory\n" );
    ourFree( buffer );
    fclose( fp );
    *status = -1;
    return( NULL );
  }

  /*
  ** Initialize values.
  */
  chromatData->fileType           = ABIFormat;
  chromatData->primerLoc          = 0;
  chromatData->avgSpacing         = 0.0;
  chromatData->machineName[0]     = '\0';
  chromatData->sampleName[0]      = '\0';
  chromatData->primerID[0]        = '\0';
  chromatData->signalStrength[0]  = 0;
  chromatData->signalStrength[1]  = 0;
  chromatData->signalStrength[2]  = 0;
  chromatData->signalStrength[3]  = 0;
  chromatData->gelName[0]         = '\0';
  chromatData->laneNumber         = 0;
  chromatData->processing[0]      = '\0';
  chromatData->reTracker[0]       = '\0';
  chromatData->comment[0]         = '\0';
  chromatData->convProg[0]        = '\0';
  pstrcpy( chromatData->fileName, fn, PHRED_PATH_MAX );
  strcpy( chromatData->source, "ABI 373A or 377" );
  memset( chromatData->thumbPrint, 0, 10 * sizeof( char ) );

  /*
  ** Read thumb print.
  */
  offset = (size_t)dataIndex[THMPRT].offset;
  numbyt = (size_t)dataIndex[THMPRT].numbyt;
  occur  =         dataIndex[THMPRT].occur;
  if( occur == -1 ||
      fseek( fp, offset, SEEK_SET ) != 0 ||
      fread( buffer, (size_t)1, numbyt, fp )   != numbyt )
  {
/*
    fprintf( stderr, "unable to read thumbprint\n" );
*/
  }
  else
  {
    memcpy( chromatData->thumbPrint, buffer, 10 );
  }

  /*
  ** Read machine name.
  */
  offset = (size_t)dataIndex[MCHNAM].offset;
  numbyt = (size_t)dataIndex[MCHNAM].numbyt;
  occur  =         dataIndex[MCHNAM].occur;
  if( occur == 1 )
  {
    /*
    ** Truncate excessively long strings.
    */
    if( numbyt > PHRED_MAX_STRING_LEN )
    {
      numbyt = PHRED_MAX_STRING_LEN;
    }

    if( numbyt > 4 )
    {
      if( fseek( fp, offset, SEEK_SET ) == 0 &&
          fread( buffer, (size_t)1, numbyt, fp ) == numbyt )
      {
        strncpy( chromatData->machineName, buffer + 1, numbyt - 1 );
        chromatData->machineName[numbyt-1] = '\0';
      }
    }
    else
    {
      strncpy( chromatData->machineName,
               (char *) &(dataIndex[MCHNAM].offset) + 1,
               numbyt - 1 );
      chromatData->machineName[numbyt-1] = '\0';
    }
  }

  /*
  ** Read dye primer.
  */
  offset = (size_t)dataIndex[DYEPRI].offset;
  numbyt = (size_t)dataIndex[DYEPRI].numbyt;
  occur  =         dataIndex[DYEPRI].occur;
  if( occur == 1 )
  {
    /*
    ** Truncate excessively long strings.
    */
    if( numbyt > PHRED_MAX_STRING_LEN )
    {
      numbyt = PHRED_MAX_STRING_LEN;
    }

    if( numbyt > 4 )
    {
      if( fseek( fp, offset, SEEK_SET ) == 0 &&
          fread( buffer, (size_t)1, numbyt, fp ) == numbyt )
      {
        strncpy( chromatData->primerID, buffer + 1, numbyt - 1 );
        chromatData->primerID[numbyt-1] = '\0';
      }
    }
    else
    {
      strncpy( chromatData->primerID,
               (char *) &(dataIndex[DYEPRI].offset) + 1,
               numbyt - 1 );
      chromatData->primerID[numbyt-1] = '\0';
    }
  }

  /*
  ** Read sample name.
  */
  offset = (size_t)dataIndex[SMPNAM].offset;
  numbyt = (size_t)dataIndex[SMPNAM].numbyt;
  occur  =         dataIndex[SMPNAM].occur;
  if( occur == 1 )
  {
    /*
    ** Truncate excessively long strings.
    */
    if( numbyt > PHRED_MAX_STRING_LEN )
    {
      numbyt = PHRED_MAX_STRING_LEN;
    }

    if( numbyt > 4 )
    {
      if( fseek( fp, offset, SEEK_SET ) == 0 &&
          fread( buffer, (size_t)1, numbyt, fp ) == numbyt )
      {
        strncpy( chromatData->sampleName, buffer + 1, numbyt - 1 );
        chromatData->sampleName[numbyt-1] = '\0';
      }
    }
    else
    {
      strncpy( chromatData->sampleName,
               (char *) &(dataIndex[SMPNAM].offset) + 1,
               numbyt - 1 );
      chromatData->sampleName[numbyt-1] = '\0';
    }
  }

  /*
  ** Read lane number
  */
  occur  =         dataIndex[LANENM].occur;
  if( occur == -1 )
  {
    fprintf( stderr, "unable to read lane number\n" );
  }
  else
  {
    chromatData->laneNumber = (int)inSwpUint2( (char *) &(dataIndex[LANENM].offset) );
  }

  /*
  ** Read gel name.
  */
  offset = (size_t)dataIndex[GELNAM].offset;
  numbyt = (size_t)dataIndex[GELNAM].numbyt;
  occur  =         dataIndex[GELNAM].occur;
  if( occur == 1 )
  {
    /*
    ** Truncate excessively long strings.
    */
    if( numbyt > PHRED_MAX_STRING_LEN )
    {
      numbyt = PHRED_MAX_STRING_LEN;
    }

    if( numbyt > 4 )
    {
      if( fseek( fp, offset, SEEK_SET ) == 0 &&
          fread( buffer, (size_t)1, numbyt, fp ) == numbyt )
      {
        strncpy( chromatData->gelName, buffer + 1, numbyt - 1 );
        chromatData->gelName[numbyt-1] = '\0';
      }
    }
    else
    {
      strncpy( chromatData->gelName,
               (char *) &(dataIndex[GELNAM].offset) + 1,
               numbyt - 1 );
      chromatData->gelName[numbyt-1] = '\0';
    }
  }

  /*
  ** Read comment.
  */
  offset = (size_t)dataIndex[COMMNT].offset;
  numbyt = (size_t)dataIndex[COMMNT].numbyt;
  occur  =         dataIndex[COMMNT].occur;
  if( occur == 1 )
  {
    /*
    ** Truncate excessively long strings.
    */
    if( numbyt > PHRED_MAX_STRING_LEN )
    {
      numbyt = PHRED_MAX_STRING_LEN;
    }

    if( numbyt > 4 )
    {
      if( fseek( fp, offset, SEEK_SET ) == 0 &&
          fread( buffer, (size_t)1, numbyt, fp ) == numbyt )
      {
        strncpy( chromatData->comment, buffer + 1, numbyt - 1 );
        chromatData->comment[numbyt-1] = '\0';
      }
    }
    else
    {
      strncpy( chromatData->comment,
               (char *) &(dataIndex[COMMNT].offset) + 1,
               numbyt - 1 );
      chromatData->comment[numbyt-1] = '\0';
    }
  }

  /*
  ** Return if there are no trace data.
  */
  if( numPoint == 0 )
  {
    fclose( fp );
    ourFree( buffer );
    *status = 2;
    return( chromatData );
  }

  /*
  ** Return if there is no map.
  */
  if( dataIndex[BASMAP].occur == -1 )
  {
    if( chromatData->numPoint > 0 )
    {
      for( i = 0; i < 4; ++i )
      {
        ourFree( (char *)chromatData->trace[i] );
        chromatData->trace[i] = NULL;
      }
    }
    chromatData->numPoint = 0;
    chromatData->numBase  = 0;
    fclose( fp );
    ourFree( buffer );
    fprintf( stderr, "readABI: no base map\n" );
    *status = 1;
    return( chromatData );
  }

  /*
  ** Get base mapping.  Each element is a letter from the
  ** set {ACGT}.
  */
  if( dataIndex[BASMAP].numbyt <= 4 )
  {
    dataIndex[BASMAP].offset = inSwpUint4( (char *) &(dataIndex[BASMAP].offset) );
    baseMap[0] = ( dataIndex[BASMAP].offset >> 24 ) & 0xff;
    baseMap[1] = ( dataIndex[BASMAP].offset >> 16 ) & 0xff;
    baseMap[2] = ( dataIndex[BASMAP].offset >>  8 ) & 0xff;
    baseMap[3] =   dataIndex[BASMAP].offset         & 0xff;
  }
  else
  {
    offset = (size_t)dataIndex[BASMAP].offset;
    numbyt = (size_t)dataIndex[BASMAP].numbyt;
    numwrd =         dataIndex[SIGSTR].numwrd;
    occur  =         dataIndex[BASMAP].occur;
    if( occur == -1 ||
        fseek( fp, offset, SEEK_SET ) != 0 ||
        fread( buffer, (size_t)1, numbyt, fp )   != numbyt )
    {
      fprintf( stderr, "readABI: error: unable to read base map\n" );
      chromatData->numPoint = 0;
      chromatData->numBase  = 0;
      fclose( fp );
      ourFree( buffer );
      *status = 1;
      return( chromatData );
    }
    else
    {
      baseMap[0] = buffer[0];
      baseMap[1] = buffer[1];
      baseMap[2] = buffer[2];
      baseMap[3] = buffer[3];
    }
  }

  /*
  ** Read signal strength.
  */
  offset = (size_t)dataIndex[SIGSTR].offset;
  numbyt = (size_t)dataIndex[SIGSTR].numbyt;
  numwrd =         dataIndex[SIGSTR].numwrd;
  occur  =         dataIndex[SIGSTR].occur;
  if( occur == -1 ||
      fseek( fp, offset, SEEK_SET ) != 0 ||
      fread( buffer, (size_t)1, numbyt, fp )   != numbyt )
  {
/*
    fprintf( stderr, "unable to read signal strength\n" );
*/
  }
  else
  {
    uint2 sig1, sig2, sig3, sig4;
    cptr = buffer;
    sig1 = inSwpUint2( cptr );
    cptr += 2;
    sig2 = inSwpUint2( cptr );
    cptr += 2;
    sig3 = inSwpUint2( cptr );
    cptr += 2;
    sig4 = inSwpUint2( cptr );
    chromatData->signalStrength[baseChar[(int)baseMap[0]]] = (int)sig1;
    chromatData->signalStrength[baseChar[(int)baseMap[1]]] = (int)sig2;
    chromatData->signalStrength[baseChar[(int)baseMap[2]]] = (int)sig3;
    chromatData->signalStrength[baseChar[(int)baseMap[3]]] = (int)sig4;
  }

  /*
  ** Read average peak spacing.
  */
  offset = dataIndex[AVGSPC].offset;
  occur  = dataIndex[AVGSPC].occur;
  if( occur == -1 )
  {
/*
    fprintf( stderr, "unable to read peak spacing\n" );
*/
  }
  else
  {
    offset = inSwpUint4( ( char *) &offset );
    chromatData->avgSpacing = int2float( offset );
  }

  /*
  ** Read primer position.
  */
  offset = dataIndex[PRIPOS].offset;
  occur  = dataIndex[PRIPOS].occur;
  if( occur == -1 )
  {
/*
    fprintf( stderr, "unable to read primer position\n" );
*/
  }
  else
  {
    offset = inSwpUint4( (char *) &offset );
    chromatData->primerLoc = (int)( offset >> 16 );
  }

  /*
  ** Read first trace.
  */
  offset = (size_t)dataIndex[TRACE1].offset;
  numbyt = (size_t)dataIndex[TRACE1].numbyt;
  numwrd =         dataIndex[TRACE1].numwrd;
  occur  =         dataIndex[TRACE1].occur;
  if( occur == -1 )
  {
    fprintf( stderr, "unable to read %c trace\n", baseMap[0] );
    fclose( fp );
    ourFree( buffer );
    chromatData->numPoint = 0;
    chromatData->numBase  = 0;
    *status = 2;
    return( chromatData );
  }
  if( numbyt > 4 )
  {
    if( fseek( fp, offset, SEEK_SET ) != 0 ||
        fread( buffer, (size_t)1, numbyt, fp )   != numbyt )
    {
      fprintf( stderr, "unable to read %c trace\n", baseMap[0] );
      fclose( fp );
      ourFree( buffer );
      chromatData->numPoint = 0;
      chromatData->numBase  = 0;
      *status = 2;
      return( chromatData );
    }
    cptr = buffer;
  }
  else
  {
    cptr = (char *) &offset;
  }
  fptr = chromatData->trace[baseChar[(int)baseMap[0]]];
  for( i = 0; i < (int)numwrd; ++i )
  {
    if( cptr[0] == '\015' || cptr[1] == '\015' )
    {
      crFlg = 1;
    }
    *fptr = (FLOAT)inSwpUint2( cptr );
    cptr += 2;
    ++fptr;
  }

  /*
  ** Read second trace.
  */
  offset = (size_t)dataIndex[TRACE2].offset;
  numbyt = (size_t)dataIndex[TRACE2].numbyt;
  numwrd =         dataIndex[TRACE2].numwrd;
  occur  =         dataIndex[TRACE2].occur;
  if( occur == -1 )
  {
    fprintf( stderr, "unable to read %c trace\n", baseMap[1] );
    fclose( fp );
    ourFree( buffer );
    chromatData->numPoint = 0;
    chromatData->numBase  = 0;
    *status = 2;
    return( chromatData );
  }
  if( numbyt > 4 )
  {
    if( fseek( fp, offset, SEEK_SET ) != 0 ||
        fread( buffer, (size_t)1, numbyt, fp )   != numbyt )
    {
      fprintf( stderr, "unable to read %c trace\n", baseMap[1] );
      fclose( fp );
      ourFree( buffer );
      chromatData->numPoint = 0;
      chromatData->numBase  = 0;
      *status = 2;
      return( chromatData );
    }
    cptr = buffer;
  }
  else
  {
    cptr = (char *) &offset;
  }
  fptr = chromatData->trace[baseChar[(int)baseMap[1]]];
  for( i = 0; i < (int)numwrd; ++i )
  {
    if( cptr[0] == '\015' || cptr[1] == '\015' )
    {
      crFlg = 1;
    }
    *fptr = (FLOAT)inSwpUint2( cptr );
    cptr += 2;
    ++fptr;
  }


  /*
  ** Read third trace.
  */
  offset = (size_t)dataIndex[TRACE3].offset;
  numbyt = (size_t)dataIndex[TRACE3].numbyt;
  numwrd =         dataIndex[TRACE3].numwrd;
  occur  =         dataIndex[TRACE3].occur;
  if( occur == -1 )
  {
    fprintf( stderr, "unable to read %c trace\n", baseMap[2] );
    fclose( fp );
    ourFree( buffer );
    chromatData->numPoint = 0;
    chromatData->numBase  = 0;
    *status = 2;
    return( chromatData );
  }
  if( numbyt > 4 )
  {
    if( fseek( fp, offset, SEEK_SET ) != 0 ||
        fread( buffer, (size_t)1, numbyt, fp )   != numbyt )
    {
      fprintf( stderr, "unable to read %c trace\n", baseMap[2] );
      fclose( fp );
      ourFree( buffer );
      chromatData->numPoint = 0;
      chromatData->numBase  = 0;
      *status = 2;
      return( chromatData );
    }
    cptr = buffer;
  }
  else
  {
    cptr = (char *) &offset;
  }
  fptr = chromatData->trace[baseChar[(int)baseMap[2]]];
  for( i = 0; i < (int)numwrd; ++i )
  {
    if( cptr[0] == '\015' || cptr[1] == '\015' )
    {
      crFlg = 1;
    }
    *fptr = (FLOAT)inSwpUint2( cptr );
    cptr += 2;
    ++fptr;
  }


  /*
  ** Read fourth trace.
  */
  offset = (size_t)dataIndex[TRACE4].offset;
  numbyt = (size_t)dataIndex[TRACE4].numbyt;
  numwrd =         dataIndex[TRACE4].numwrd;
  occur  =         dataIndex[TRACE4].occur;
  if( occur == -1 )
  {
    fprintf( stderr, "unable to read %c trace\n", baseMap[3] );
    fclose( fp );
    ourFree( buffer );
    chromatData->numPoint = 0;
    chromatData->numBase  = 0;
    *status = 2;
    return( chromatData );
  }
  if( numbyt > 4 )
  {
    if( fseek( fp, offset, SEEK_SET ) != 0 ||
        fread( buffer, (size_t)1, numbyt, fp )   != numbyt )
    {
      fprintf( stderr, "unable to read %c trace\n", baseMap[3] );
      fclose( fp );
      ourFree( buffer );
      chromatData->numPoint = 0;
      chromatData->numBase  = 0;
      *status = 2;
      return( chromatData );
    }
    cptr = buffer;
  }
  else
  {
    cptr = (char *) &offset;
  }
  fptr = chromatData->trace[baseChar[(int)baseMap[3]]];
  for( i = 0; i < (int)numwrd; ++i )
  {
    if( cptr[0] == '\015' || cptr[1] == '\015' )
    {
      crFlg = 1;
    }
    *fptr = (FLOAT)inSwpUint2( cptr );
    cptr += 2;
    ++fptr;
  }

  /*
  ** Find minimum and maximum trace values.
  */
  findTraceExtrema( chromatData );

  /*
  ** Read bases.
  */
  offset = (size_t)dataIndex[BASES].offset;
  numbyt = (size_t)dataIndex[BASES].numbyt;
  numwrd =         dataIndex[BASES].numwrd;
  occur  =         dataIndex[BASES].occur;
  if( occur == -1 )
  {
    fprintf( stderr, "unable to read bases\n" );
    fclose( fp );
    ourFree( buffer );
    chromatData->numBase  = 0;
    *status = 3;
    return( chromatData );
  }
  if( numbyt > 4 )
  {
    if( fseek( fp, offset, SEEK_SET ) != 0 ||
        fread( buffer, (size_t)1, numbyt, fp )   != numbyt )
    {
      fprintf( stderr, "unable to read bases\n" );
      fclose( fp );
      ourFree( buffer );
      chromatData->numBase  = 0;
      *status = 3;
      return( chromatData );
    }
    cptr = buffer;
  }
  else
  {
    cptr = (char *)&offset;
  }
  for( i = 0; i < (int)numwrd; ++i )
  {
    chromatData->base[i]     = *cptr;
    chromatData->baseQual[i] = 0;
    ++cptr;
  }

  /*
  ** Read base positions.
  */
  offset = (size_t)dataIndex[BASPOS].offset;
  numbyt = (size_t)dataIndex[BASPOS].numbyt;
  numwrd =         dataIndex[BASPOS].numwrd;
  occur  =         dataIndex[BASPOS].occur;
  if( occur == -1 )
  {
    fprintf( stderr, "unable to read base positions\n" );
    fclose( fp );
    ourFree( buffer );
    chromatData->numBase  = 0;
    *status = 3;
    return( chromatData );
  }
  if( numbyt > 4 )
  {
    if( fseek( fp, offset, SEEK_SET ) != 0 ||
        fread( buffer, (size_t)1, numbyt, fp )   != numbyt )
    {
      fprintf( stderr, "unable to read base positions\n" );
      fclose( fp );
      ourFree( buffer );
      chromatData->numBase  = 0;
      *status = 3;
      return( chromatData );
    }
    cptr = buffer;
  }
  else
  {
    cptr = (char *) &offset;
  }
  for( i = 0; i < (int)numwrd; ++i )
  {
    chromatData->baseLoc[i]  = (int)inSwpUint2( cptr );
    chromatData->baseQual[i] = 0;
    cptr += 2;
  }

  /*
  ** Free memory.
  */
  ourFree( buffer );

  /*
  ** Finished with the file
  */
  fclose( fp );

  /*
  ** Possibly transferred data in ASCII mode if <cr>s absent.
  */
  if( crFlg == 0 &&
      option->verboseOption == 1 && option->verboseLevel >= 1 )
  {
    fprintf( stderr,
             "readABI: corrupted file? - FTP in binary mode\n" );
  }

  if( option->verboseOption && option->verboseLevel >= 16 )
  {
    fprintf( stderr,
             "readABI: machine name:     %s\n",
             chromatData->machineName );
    fprintf( stderr,
             "readABI: gel name:         %s\n",
             chromatData->gelName );
    fprintf( stderr,
             "readABI: sample name:      %s\n",
             chromatData->sampleName );
    fprintf( stderr,
             "readABI: primer ID:        %s\n",
             chromatData->primerID );
    fprintf( stderr,
             "readABI: lane number:      %d\n",
             chromatData->laneNumber );
    fprintf( stderr,
             "readABI: comment:          %s\n",
             chromatData->comment );
    fprintf( stderr,
             "readABI: number of scans:  %d\n",
             chromatData->numPoint );
    fprintf( stderr,
             "readABI: base order:       %c %c %c %c\n",
             baseMap[0],
             baseMap[1],
             baseMap[2],
             baseMap[3] );
  }

  *status = 0;

  return( chromatData );
}
