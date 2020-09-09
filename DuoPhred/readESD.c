/** readESD.c **/

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
** Read Molecular Dynamics ESD files.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <sys/types.h>
#include <sys/stat.h>

#include "phred.h"
#include "readESD.h"


#define TRAILER_SIZE	632
#define POINT_SIZE	20

/*
** Intel machines write ESD files so we
** must deal with little-endian words.
*/
#ifdef ANSI_C
static int4 inSwpSint4L( unsigned char *ptr )
#else
static int4 inSwpSint4L( ptr )
unsigned char *ptr;
#endif
{
  int4 i;

  ptr += 3;

  i = *ptr & 0xff;
  i = i << 8;
  --ptr;

  i |= *ptr & 0xff;
  i = i << 8;
  --ptr;

  i |= *ptr & 0xff;
  i = i << 8;
  --ptr;

  i |= *ptr & 0xff;

  return( i );
}

#ifdef ANSI_C
static uint2 inSwpUint2L( unsigned char *ptr )
#else
static uint2 inSwpUint2L( ptr )
unsigned char *ptr;
#endif
{
  uint2 i;

  ptr += 1;

  i = *ptr & 0xff;
  i = i << 8;
  --ptr;

  i |= *ptr & 0xff;

  return( i );
}

#ifdef ANSI_C
static uint4 inSwpUint3L( unsigned char *ptr )
#else
static uint4 inSwpUint3L( ptr )
unsigned char *ptr;
#endif
{
  uint4 i;

  ptr += 2;

  i = *ptr & 0xff;
  i = i << 8;
  --ptr;

  i |= *ptr & 0xff;
  i = i << 8;
  --ptr;

  i |= *ptr & 0xff;

  return( i );
}

#ifdef ANSI_C
static uint4 inSwpUint4L( unsigned char *ptr )
#else
static uint4 inSwpUint4L( ptr )
unsigned char *ptr;
#endif
{
  uint4 i;

  ptr += 3;

  i = *ptr & 0xff;
  i = i << 8;
  --ptr;

  i |= *ptr & 0xff;
  i = i << 8;
  --ptr;

  i |= *ptr & 0xff;
  i = i << 8;
  --ptr;

  i |= *ptr & 0xff;

  return( i );
}

#ifdef ANSI_C
static float inSwpFloatL( unsigned char *ptr )
#else
static float inSwpFloatL( ptr )
unsigned char *ptr;
#endif
{
  unsigned int r;

  ptr += 3;

  r = *ptr & 0xff;
  r = r << 8;
  --ptr;

  r |= *ptr & 0xff;
  r = r << 8;
  --ptr;

  r |= *ptr & 0xff;
  r = r << 8;
  --ptr;

  r |= *ptr & 0xff;

  return( *(float *) &r );
}


/*
** Extract information from trailer block.
*/
#ifdef ANSI_C
static int extractTrailerData( unsigned char *buf, ESDTrailer *trailer )
#else
static int extractTrailerData( buf, trailer )
unsigned char *buf;
ESDTrailer *trailer;
#endif
{
  unsigned char *cptr;

  cptr = buf;

  strncpy( (char *)trailer->m_tcSampleName, (char *)cptr, 24 );
  cptr += 24;

  strncpy( (char *)trailer->m_tcOperator, (char *)cptr, 30 );
  cptr += 30;

  trailer->m_fRunVoltage = inSwpFloatL( cptr );
  cptr += 4;

  trailer->m_fRunTime = inSwpFloatL( cptr );
  cptr += 4;

  trailer->m_fInjectionTime = inSwpFloatL( cptr );
  cptr += 4;

  trailer->m_fInjectionVoltage = inSwpFloatL( cptr );
  cptr += 4;

  trailer->m_fTemperature = inSwpFloatL( cptr );
  cptr += 4;

  trailer->m_nLazerIndex = inSwpSint4L( cptr );
  cptr += 4;

  trailer->m_fPMTVoltage1 = inSwpFloatL( cptr );
  cptr += 4;

  trailer->m_fPMTVoltage2 = inSwpFloatL( cptr );
  cptr += 4;

  strncpy( (char *)trailer->m_csBeamSplitterA, (char *)cptr, 30 );
  cptr += 30;

  strncpy( (char *)trailer->m_csBeamSplitterB, (char *)cptr, 30 );
  cptr += 30;

  strncpy( (char *)trailer->m_csFilter1, (char *)cptr, 30 );
  cptr += 30;
  
  strncpy( (char *)trailer->m_csFilter2, (char *)cptr, 30 );
  cptr += 30;
 
  strncpy( (char *)trailer->m_csFilter3, (char *)cptr, 30 );
  cptr += 30;
 
  strncpy( (char *)trailer->m_csFilter4, (char *)cptr, 30 );
  cptr += 30;
 
  trailer->m_fFillTime = inSwpFloatL( cptr );
  cptr += 4;

  trailer->m_fFlushTime1 = inSwpFloatL( cptr );
  cptr += 4;

  trailer->m_fFlushTime2 = inSwpFloatL( cptr );
  cptr += 4;

  trailer->m_fFlushTime3 = inSwpFloatL( cptr );
  cptr += 4;

  trailer->m_fPrerunVoltage = inSwpFloatL( cptr );
  cptr += 4;

  trailer->m_fRelaxationTime = inSwpFloatL( cptr );
  cptr += 4;

  trailer->m_fPrerunTime = inSwpFloatL( cptr );
  cptr += 4;

  cptr += 64;

  strncpy( (char *)trailer->m_csBaseOrder[0], (char *)cptr, 10 );
  cptr += 10;

  strncpy( (char *)trailer->m_csBaseOrder[1], (char *)cptr, 10 );
  cptr += 10;

  strncpy( (char *)trailer->m_csBaseOrder[2], (char *)cptr, 10 );
  cptr += 10;

  strncpy( (char *)trailer->m_csBaseOrder[3], (char *)cptr, 10 );
  cptr += 10;

  strncpy( (char *)trailer->m_InstID, (char *)cptr, 22 );
  cptr += 22;

  strncpy( (char *)trailer->m_csParamSetName, (char *)cptr, 30 );
  cptr += 30;

  trailer->m_dMagicNumber = inSwpSint4L( cptr );
  cptr += 4;

  trailer->m_dExtHeaderOffset = inSwpSint4L( cptr );
  cptr += 4;

  trailer->m_dExtHeaderLen = inSwpSint4L( cptr );
  cptr += 4;

  trailer->m_fPreinjectionTime = inSwpFloatL( cptr );
  cptr += 4;

  trailer->m_fPreinjectionVoltage = inSwpFloatL( cptr );
  cptr += 4;

  trailer->m_fLowPressureTime = inSwpFloatL( cptr );
  cptr += 4;

  trailer->m_fUserInputTime = inSwpFloatL( cptr );
  cptr += 4;

  trailer->m_fSleepTime = inSwpFloatL( cptr );
  cptr += 4;

  trailer->m_fSleepTemperature = inSwpFloatL( cptr );
  cptr += 4;

  strncpy( (char *)trailer->m_byReserved, (char *)cptr, 82 );
  cptr += 82;

  strncpy( (char *)trailer->m_csRunID, (char *)cptr, 6 );
  cptr += 6;

  strncpy( (char *)trailer->m_csPlateID, (char *)cptr, 18 );
  cptr += 18;

  trailer->m_dwPrevRunID = inSwpSint4L( cptr );
  cptr += 4;

  trailer->m_bDWORDFile = inSwpSint4L( cptr );
  cptr += 4;

  trailer->m_bSaturated = inSwpSint4L( cptr );
  cptr += 4;

  trailer->m_fVersion = inSwpFloatL( cptr );
  cptr += 4;

  trailer->m_tRunStartTime = inSwpSint4L( cptr );
  cptr += 4;

  trailer->m_tRunStopTime = inSwpSint4L( cptr );
  cptr += 4;

  trailer->m_dwNumberOfLines = inSwpSint4L( cptr );
  cptr += 4;

  trailer->m_nScanRate = inSwpSint4L( cptr );
  cptr += 4;

  trailer->m_byFileType = *cptr;
  cptr += 1;

  trailer->m_byNumberOfCapillaries = *cptr;
  cptr += 1;

  trailer->m_byCapUsage = *cptr;
  cptr += 1;

  trailer->m_byDummy = *cptr;
  cptr += 1;

  trailer->m_uiSecurityCode = inSwpUint4L( cptr );

  return( 0 );

}

#define LABEL_ELEMENT	0
#define VALUE_ELEMENT	1



/*
** Process extended header lists recursively.
*/
#ifdef ANSI_C
int recurList( unsigned char *buf, long int lenData, int numElement,
               long int *lib, char *inLabel,
               ESDExtHeader *header )
#else
int recurList( buf, lenData, numElement, lib, inLabel, header )
unsigned char *buf;
long int lenData;
int numElement;
long int *lib;
char *inLabel;
ESDExtHeader *header;
#endif
{
  int ie, k;
  long int ib;
  int elementType;
  int format;
  int lenSize;
  int lenByte;
  unsigned char headerByte;
  char label[3*ESD_STR_MAX];

  ib = *lib;

  ie = 1;
  elementType = LABEL_ELEMENT;
  while( ib < lenData && ie <= numElement )
  {
    headerByte = buf[ib];
    ++ib;

    format = headerByte >> 2;
    lenSize = headerByte & 0x03;

    /*
    ** Watch for end of extended header data block.
    */
    if( lenSize == 0 && format == 077 )
    {
      break;
    }

    if( lenSize == 1 )
    {
      lenByte = buf[ib] & 0xff;
      ++ib;
    }
    else
    if( lenSize == 2 )
    {
      lenByte = inSwpUint2L( &(buf[ib]) );
      ib += 2;
    }
    else
    if( lenSize == 3 )
    {
      lenByte = inSwpUint3L( &(buf[ib]) );
      ib += 3;
    }
    else
    {
      fprintf( stderr,
               "readESD: unsupported header length\n" );
      return( ERROR );
    }

    /*
    ** Are we expecting a label?
    ** The value should be a string.
    */
    if( elementType == LABEL_ELEMENT )
    {
/*
      if( format == 01 || format == 02 )
*/
      if( format == 01 )
      {
        strcpy( label, inLabel );
        strcat( label, "." );
        strncat( label, (char *) &(buf[ib]), lenByte );
      }
      else
      if( format == 02 )
      {
      }
      else
      {
        fprintf( stderr,
                 "readESD: recurs: %s: unable to extract extended header: bad format: %d\n",
                 inLabel, format );
        *lib = ib;
        return( ERROR );
      }
      elementType = VALUE_ELEMENT;
      ib += lenByte;
    }
    else
    {
/*
      if( format == 01 || format == 02 )
*/
      if( format == 01 )
      {
        /*
        ** copy string.
        */
        if( strcmp( label, "CHEMISTRY.NAME" ) == 0 )
        {
          k = lenByte < ESD_STR_MAX ? lenByte : ESD_STR_MAX - 1;
          strncpy( header->primerName, (char *) &(buf[ib]), k );
        }
        else
        if( strcmp( label, "CHEMISTRY.CHANNEL1.BASE" ) == 0 )
        {
          header->baseMap[0] = buf[ib];
        }
        else
        if( strcmp( label, "CHEMISTRY.CHANNEL2.BASE" ) == 0 )
        {
          header->baseMap[1] = buf[ib];
        }
        else
        if( strcmp( label, "CHEMISTRY.CHANNEL3.BASE" ) == 0 )
        {
          header->baseMap[2] = buf[ib];
        }
        else
        if( strcmp( label, "CHEMISTRY.CHANNEL4.BASE" ) == 0 )
        {
          header->baseMap[3] = buf[ib];
        }

        elementType = LABEL_ELEMENT;
        ib += lenByte;
      }
      else
      if( format == 00 )
      {
        /*
        ** Process list.
        */
        recurList( buf, lenData, lenByte, &ib, label, header );
        elementType = LABEL_ELEMENT;
      }
      else
      {
        /*
        ** Get the value if necessary.
        */
        elementType = LABEL_ELEMENT;
        ib += lenByte;
      }
    }

    ++ie;
  }

  *lib = ib;

  return( OK );
}



/*
** Read extended header block.
*/
#ifdef ANSI_C
int extractExtendedHeader( unsigned char *buf, long int len, ESDExtHeader *header )
#else
int extractExtendedHeader( buf, len, header )
unsigned char *buf;
long int len;
ESDExtHeader *header;
#endif
{
  int k, m;
  int format;
  int lenSize;
  int lenByte;
  int numByte;
  int elementType;
  long int ib;
  unsigned char headerByte;
  char label[ESD_STR_MAX];

  ib = 0;
  elementType = LABEL_ELEMENT;
  while( ib < len )
  {
    headerByte = buf[ib];
    ++ib;

    format = headerByte >> 2;
    lenSize = headerByte & 0x03;

    /*
    ** Watch for end of extended header data block.
    */
    if( lenSize == 0 && format == 077 )
    {
      break;
    }

    if( lenSize == 1 )
    {
      lenByte = buf[ib] & 0xff;
      ++ib;
    }
    else
    if( lenSize == 2 )
    {
      lenByte = inSwpUint2L( &(buf[ib]) );
      ib += 2;
    }
    else
    if( lenSize == 3 )
    {
      lenByte = inSwpUint3L( &(buf[ib]) );
      ib += 3;
    }
    else
    {
      fprintf( stderr,
               "readESD: unsupported header length\n" );
      return( ERROR );
    }

    /*
    ** Are we expecting a label?
    ** The value should be a string.
    */
    if( elementType == LABEL_ELEMENT )
    {
/*
      if( format == 01 || format == 02 )
*/
      if( format == 01 )
      {
        strncpy( label,
                 (char *)&(buf[ib]),
                 lenByte < ESD_STR_MAX ? lenByte : ( ESD_STR_MAX - 1 ) );
        label[ESD_STR_MAX-1] = '\0';
      }
      else
      if( format == 02 )
      {

      }
      else
      {
        fprintf( stderr,
                 "readESD: unable to extract extended header: element type: %d   format: %d\n",
                 elementType, format );
        return( ERROR );
      }
      elementType = VALUE_ELEMENT;
      ib += lenByte;
    }
    else
    {
/*
      if( format == 01 || format == 02 )
*/
      if( format == 01 )
      {
        /*
        ** Copy the strings.
        */
        if( strcmp( label, "PLATE ID" ) == 0 )
        {
          k = lenByte < ESD_STR_MAX ? lenByte : ESD_STR_MAX - 1;
          strncpy( header->plateID, (char *) &(buf[ib]), k );
          header->plateID[ESD_STR_MAX-1] = '\0';
        }
        else
        if( strcmp( label, "SAMPLE NAME" ) == 0 )
        {
          k = lenByte < ESD_STR_MAX ? lenByte : ESD_STR_MAX - 1;
          strncpy( header->sampleName, (char *) &(buf[ib]), k );
          header->sampleName[ESD_STR_MAX-1] = '\0';
        }
        else
        if( strcmp( label, "WELL ID" ) == 0 )
        {
          k = lenByte < ESD_STR_MAX ? lenByte : ESD_STR_MAX - 1;
          strncpy( header->wellID, (char *) &(buf[ib]), k );
          header->wellID[ESD_STR_MAX-1] = '\0';
        }
        else
        if( strcmp( label, "COMMENT" ) == 0 )
        {
          k = lenByte < ESD_STR_MAX ? lenByte : ESD_STR_MAX - 1;
          strncpy( header->comment, (char *) &(buf[ib]), k );
          header->comment[ESD_STR_MAX-1] = '\0';
        }
        else
        if( strcmp( label, "MACHINE ID" ) == 0 )
        {
          k = lenByte < ESD_STR_MAX ? lenByte : ESD_STR_MAX - 1;
          strncpy( header->machineName, (char *) &(buf[ib]), k );
          header->machineName[ESD_STR_MAX-1] = '\0';
        }
        else
        if( strcmp( label, "SEQUENCE" ) == 0 )
        {
          header->numBase = lenByte - 1;
          header->base = (char *)malloc( lenByte );
          if( header->base == NULL )
          {
            fprintf( stderr,
                     "readESD: unable to allocate memory\n" );
            return( ERROR );
          }
          for( k = 0; k < header->numBase; ++k )
          {
            header->base[k] = buf[ib+k];
          }
        }

        elementType = LABEL_ELEMENT;
        ib += lenByte;
      }
      else
      if( format == 00 )
      {
        /*
        ** Process list.
        */
        recurList( buf, len, lenByte, &ib, label, header );
        elementType = LABEL_ELEMENT;
      }
      else
      {
        /*
        ** Get the value if necessary.
        */
        if( strcmp( label, "PEAK POSITIONS" ) == 0 )
        {
          if( format == 054 )
          {
            m = lenByte / 4;
            numByte = m * sizeof( int );
            header->baseLoc = (int *)malloc( numByte );
            if( header->baseLoc == NULL )
            {
              fprintf( stderr,
                       "readESD: unable to allocate memory\n" );
              return( ERROR );
            }
            for( k = 0; k < m; ++k )
            {
              header->baseLoc[k] = inSwpSint4L( &(buf[ib]) );
              ib += 4;
            }
          }
          else
          {
            fprintf( stderr,
                     "readESD: unexpected peak position data format\n" );
            ib += lenByte;
          }
          elementType = LABEL_ELEMENT;
        }
        else
        if( strcmp( label, "START TRACE PROCESSING" ) == 0 )
        {
          if( format == 054 )
          {
            header->baseLocOffset = inSwpSint4L( &(buf[ib]) );
          }
          else
          {
            fprintf( stderr,
                     "readESD: unexpected start trace processing data format\n" );
          }
          ib += lenByte;
          elementType = LABEL_ELEMENT;
        }
        else
        {
          elementType = LABEL_ELEMENT;
          ib += lenByte;
        }
      }
    }
  }

  return( OK );
}



/*
** Check whether a file is likely an ESD file.
*/
#ifdef ANSI_C
int isESD( char *filename )
#else
int isESD( filename )
char *filename;
#endif
{
  int istat;
  int nitem;
  long fileSize;
  struct STAT_STRUCT statBuf;
  unsigned char trailerBuf[TRAILER_SIZE];
  ESDTrailer trailer;
  FILE *fp;
  Option *option;

  option = getOption();

  /*
  ** Get file size.
  */
  istat = FILE_STATUS( filename, &statBuf );
  if( istat == -1 )
  {
    if( option->verboseOption && option->verboseLevel >= 32 )
    {
      fprintf( stderr,
               "isESD: unable to `stat' file: %s\n",
               filename );
    }
    return( 0 );
  }
  fileSize = statBuf.st_size;

  /*
  ** Open file.
  */
  fp = fopen( filename, "rb" );
  if( fp == NULL )
  {
    if( option->verboseOption && option->verboseLevel >= 32 )
    {
      fprintf( stderr,
               "isESD: unable to open file: %s\n",
               filename );
    }
    return( 0 );
  }

  /*
  ** Move to trailer.
  */
  if( fseek( fp, (long) -TRAILER_SIZE, 2 ) != 0 )
  {
    if( option->verboseOption && option->verboseLevel >= 32 )
    {
      fprintf( stderr,
               "isESD: unable to `seek' file trailer: %s\n",
               filename );
    }
    fclose( fp );
    return( 0 );
  }

  /*
  ** Read trailer.
  */
  nitem = fread( trailerBuf, TRAILER_SIZE, 1, fp );
  if( nitem != 1 )
  {
    if( option->verboseOption && option->verboseLevel >= 32 )
    {
      fprintf( stderr,
               "isESD: unable to read file trailer: %s\n",
               filename );
    }
    fclose( fp );
    return( 0 );
  }

  /*
  ** Extract trailer data.
  */
  istat = extractTrailerData( trailerBuf, &trailer );
  if( istat != 0 )
  {
    if( option->verboseOption && option->verboseLevel >= 32 )
    {
      fprintf( stderr,
               "isESD: unable to extract trailer data: %s\n",
               filename );
    }
    fclose( fp );
    return( 0 );
  }

  /*
  ** Is it a raw ESD file?
  */
  if( trailer.m_bDWORDFile == 1 )
  {
    if( option->verboseOption && option->verboseLevel >= 32 )
    {
      fprintf( stderr,
               "isESD: raw ESD fle: %s\n",
               filename );
    }
    fclose( fp );
    return( 0 );
  }

  /*
  ** Check that the number of scans makes any sense.
  */
  if( fileSize - trailer.m_dwNumberOfLines * sizeof( ProcessedScan ) ==
      TRAILER_SIZE )
  {
    if( option->verboseOption && option->verboseLevel >= 32 )
    {
      fprintf( stderr,
               "isESD: ESD format without extended trailer: %s\n",
               filename );
    }
    fclose( fp );
    return( 1 );
  }
  else
  if( ( fileSize - trailer.m_dwNumberOfLines * sizeof( ProcessedScan ) ==
        trailer.m_dExtHeaderOffset ) &&
      trailer.m_dMagicNumber == 826561349 )
  {
    if( option->verboseOption && option->verboseLevel >= 32 )
    {
      fprintf( stderr,
               "isESD: ESD format with extended trailer: %s\n",
               filename );
    }
    fclose( fp );
    return( 2 );
  }

  if( option->verboseOption && option->verboseLevel >= 32 )
  {
    fprintf( stderr,
             "isESD: file size and scan number inconsistent: %s\n",
             filename );
  }

  fclose( fp );

  return( 0 );
}



/*
** status:
**
**   0 = OK
**   1 = file reading error
**   2 = no trace (and no bases assumed)
**   3 = no bases (but there is trace)
**  -1 = fatal error
*/

#ifdef ANSI_C
ChromatData *readESD( char *filename, int version, int *status )
#else
ChromatData *readESD( filename, version, status )
char *filename;
int version;
int *status;
#endif
{
  int i;
  int numByte;
  int nitem;
  int istat;
  int numPoint;
  int numBase;
  int baseChar[128];
  unsigned char *cptr;
  unsigned char *buf;
  unsigned char trailerBuf[TRAILER_SIZE];
  unsigned char pbuf[POINT_SIZE];
  char baseMap[4];
  long int extHeaderOffset;
  long int extHeaderLen;
  ChromatData *chromatData;
  ESDTrailer trailer;
  ESDExtHeader extHeader;
  FILE *fp;
  Option *option;

  option = getOption();

  /*
  ** Initialize some arrays.
  */
  baseChar[(int)'A'] = 0;
  baseChar[(int)'C'] = 1;
  baseChar[(int)'G'] = 2;
  baseChar[(int)'T'] = 3;

  baseMap[0] = (char)0;
  baseMap[1] = (char)0;
  baseMap[2] = (char)0;
  baseMap[3] = (char)0;

  /*
  ** Open file.
  */
  fp = fopen( filename, "rb" );
  if( fp == NULL )
  {
    fprintf( stderr,
             "readESD: unable to open file %s\n",
             filename );
    *status = 1;
    return( NULL );
  }

  istat = fseek( fp, (long) -TRAILER_SIZE, 2 );
  if( istat != 0 )
  {
    fprintf( stderr,
             "readESD: unable to find trailer in %s\n", filename );
    numBase = 0;
    numPoint = 0;
    chromatData = allocChromatData( numPoint, numBase );
    if( chromatData == NULL )
    {
      fprintf( stderr, "readESD: unable to allocate memory\n" );
      fclose( fp );
      *status = -1;
      return( NULL );
    }
    chromatData->fileType           = version;
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
    pstrcpy( chromatData->fileName, filename, PHRED_PATH_MAX );
    strcpy( chromatData->source, "MegaBACE" );
    memset( chromatData->thumbPrint, 0, 10 * sizeof( char ) );

    fclose( fp );
    *status = 1;
    return( chromatData );
  }

  /*
  ** Read trailer.
  */
  nitem = fread( trailerBuf, TRAILER_SIZE, 1, fp );
  if( nitem != 1 )
  {
    fprintf( stderr,
             "readESD: unable to read trailer in file %s\n", filename );
    numBase = 0;
    numPoint = 0;
    chromatData = allocChromatData( numPoint, numBase );
    if( chromatData == NULL )
    {
      fprintf( stderr, "readESD: unable to allocate memory\n" );
      fclose( fp );
      *status = -1;
      return( NULL );
    }
    chromatData->fileType           = version;
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
    pstrcpy( chromatData->fileName, filename, PHRED_PATH_MAX );
    strcpy( chromatData->source, "MegaBACE" );
    memset( chromatData->thumbPrint, 0, 10 * sizeof( char ) );

    fclose( fp );
    *status = 1;
    return( chromatData );
  }

  /*
  ** Extract trailer data.
  */
  istat = extractTrailerData( trailerBuf, &trailer );
  if( istat != 0 )
  {
    fprintf( stderr,
             "readESD: bad status: extractTrailerData in file %s\n", filename );
    numBase = 0;
    numPoint = 0;
    chromatData = allocChromatData( numPoint, numBase );
    if( chromatData == NULL )
    {
      fprintf( stderr, "readESD: unable to allocate memory\n" );
      fclose( fp );
      *status = -1;
      return( NULL );
    }
    chromatData->fileType           = version;
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
    pstrcpy( chromatData->fileName, filename, PHRED_PATH_MAX );
    strcpy( chromatData->source, "MegaBACE" );
    memset( chromatData->thumbPrint, 0, 10 * sizeof( char ) );

    fclose( fp );
    *status = 1;
    return( chromatData );
  }

  /*
  ** Copy base map.
  */
  baseMap[0] = trailer.m_csBaseOrder[0][0];
  baseMap[1] = trailer.m_csBaseOrder[1][0];
  baseMap[2] = trailer.m_csBaseOrder[2][0];
  baseMap[3] = trailer.m_csBaseOrder[3][0];

  if( option->verboseOption && option->verboseLevel >= 32 )
  {
    fprintf( stderr,
             "readESD: standard trailer: ESD format version: %f\n",
             trailer.m_fVersion );
    fprintf( stderr,
             "readESD: standard trailer: base order: %c %c %c %c\n",
             baseMap[0],
             baseMap[1],
             baseMap[2],
             baseMap[3] );
    fprintf( stderr,
             "readESD: standard trailer: number of scans: %d\n",
             trailer.m_dwNumberOfLines );
  }

  /*
  ** Number of points (scans) in file.
  */
  numPoint = trailer.m_dwNumberOfLines;

  /*
  ** Allocate chromatData structure array.
  ** (no called bases in ESD file)
  */
  numBase = 0;
  chromatData = allocChromatData( numPoint, numBase );
  if( chromatData == NULL )
  {
    fprintf( stderr, "readESD: unable to allocate memory\n" );
    fclose( fp );
    *status = -1;
    return( NULL );
  }
  chromatData->fileType           = version;
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
  pstrcpy( chromatData->fileName, filename, PHRED_PATH_MAX );
  strcpy( chromatData->source, "MegaBACE" );
  memset( chromatData->thumbPrint, 0, 10 * sizeof( char ) );

  /*
  ** Copy useful trailer data into chromatData.
  */
  pstrcpy( chromatData->sampleName, (char *)trailer.m_tcSampleName, PHRED_MAX_STRING_LEN );
  pstrcpy( chromatData->gelName,    (char *)trailer.m_csRunID,      PHRED_MAX_STRING_LEN );

  /*
  ** Read extended header if version 2.00 file.
  */
  if( version == MD2Format )
  {
    extHeader.sampleName[0] = '\0';
    extHeader.plateID[0]    = '\0';
    extHeader.primerName[0] = '\0';
    extHeader.comment[0]    = '\0';
    extHeader.numBase       = 0;

    extHeaderOffset = trailer.m_dExtHeaderOffset;
    extHeaderLen    = trailer.m_dExtHeaderLen;
    istat = fseek( fp, -extHeaderOffset, 2 );
    if( istat != 0 )
    {
      fprintf( stderr,
               "readESD: unable to find extended header\n" );
      fclose( fp );
      chromatData->numPoint = 0;
      chromatData->numBase  = 0;
      *status = 1;
      return( chromatData );
    }

    numByte = extHeaderLen * sizeof( char );
    buf = (unsigned char *)malloc( numByte );
    if( buf == NULL )
    {
      fprintf( stderr,
               "readESD: unable to allocate memory\n" );
      fclose( fp );
      chromatData->numPoint = 0;
      chromatData->numBase  = 0;
      *status = 1;
      return( chromatData );
    }

    nitem = fread( buf, 1, numByte, fp );
    if( nitem != numByte )
    {
      fprintf( stderr,
               "readESD: unable to read extended header\n" );
      fclose( fp );
      free( buf );
      chromatData->numPoint = 0;
      chromatData->numBase  = 0;
      *status = 1;
      return( chromatData );
    }

    if( extractExtendedHeader( buf, numByte, &extHeader ) == ERROR )
    {
      fclose( fp );
      free( buf );
      chromatData->numPoint = 0;
      chromatData->numBase  = 0;
      *status = -1;
      return( chromatData );
    }

    /*
    ** Copy some extended header information.
    */
    pstrcpy( chromatData->primerID,    extHeader.primerName,  PHRED_MAX_STRING_LEN );
    pstrcpy( chromatData->gelName,     extHeader.plateID,     PHRED_MAX_STRING_LEN );
    pstrcpy( chromatData->comment,     extHeader.comment,     PHRED_MAX_STRING_LEN );
    pstrcpy( chromatData->sampleName,  extHeader.sampleName,  PHRED_MAX_STRING_LEN );
    pstrcpy( chromatData->machineName, extHeader.machineName, PHRED_MAX_STRING_LEN );

    if( option->verboseOption && option->verboseLevel >= 32 )
    {
      fprintf( stderr, "readESD: extended trailer: primerID: %s\n",    extHeader.primerName );
      fprintf( stderr, "readESD: extended trailer: gelName: %s\n",     extHeader.plateID );
      fprintf( stderr, "readESD: extended trailer: comment: %s\n",     extHeader.comment );
      fprintf( stderr, "readESD: extended trailer: sampleName: %s\n",  extHeader.sampleName );
      fprintf( stderr, "readESD: extended trailer: machineName: %s\n", extHeader.machineName );
    }

    chromatData->numBase = extHeader.numBase;
    if( extHeader.numBase > 0 )
    {
      free( chromatData->base );
      chromatData->base = extHeader.base;
      free( chromatData->baseLoc );
      chromatData->baseLoc = extHeader.baseLoc;
      free( chromatData->baseQual );
      numByte = extHeader.numBase * sizeof( int );
      chromatData->baseQual = (int *)malloc( numByte );
      if( chromatData->baseQual == NULL )
      {
        fprintf( stderr,
                 "readESD: unable to allocate memory\n" );
        fclose( fp );
        chromatData->numPoint = 0;
        chromatData->numBase  = 0;
        *status = 1;
        return( chromatData );
      }
      for( i = 0; i < extHeader.numBase; ++i )
      {
        chromatData->baseLoc[i]  -= extHeader.baseLocOffset + 1;
        chromatData->baseQual[i]  = 0;
      }
    }

    if( option->verboseOption && option->verboseLevel >= 16 )
    {
      fprintf( stderr,
               "readESD: extended trailer: base order: %c %c %c %c\n",
               extHeader.baseMap[0],
               extHeader.baseMap[1],
               extHeader.baseMap[2],
               extHeader.baseMap[3] );
    }

    /*
    ** Perform some sanity checks. The base map values
    ** in the standard and extended headers must match.
    */
    if( extHeader.baseMap[0] != baseMap[0] ||
        extHeader.baseMap[1] != baseMap[1] ||
        extHeader.baseMap[2] != baseMap[2] ||
        extHeader.baseMap[3] != baseMap[3] )
    {
      fprintf( stderr,
               "readESD: inconsistent base map header information\n" );
      chromatData->numPoint = 0;
      chromatData->numBase  = 0;
      fclose( fp );
      free( buf );
      *status = 1;
      return( chromatData );
    }

    /*
    ** Free buffer memory.
    */
    free( buf );
  }

  /*
  ** Check that the baseMap has meaningful values.
  */
  if( !( baseMap[0] == 'A' || baseMap[0] == 'C' ||
         baseMap[0] == 'G' || baseMap[0] == 'T' ) ||
      !( baseMap[1] == 'A' || baseMap[1] == 'C' ||
         baseMap[1] == 'G' || baseMap[1] == 'T' ) ||
      !( baseMap[2] == 'A' || baseMap[2] == 'C' ||
         baseMap[2] == 'G' || baseMap[2] == 'T' ) ||
      !( baseMap[3] == 'A' || baseMap[3] == 'C' ||
         baseMap[3] == 'G' || baseMap[3] == 'T' )  )
  {
    fprintf( stderr,
             "readESD: bad base map value(s)\n" );
    fclose( fp );
    chromatData->numPoint = 0;
    chromatData->numBase  = 0;
    *status = 1;
    return( chromatData );
  }

  /*
  ** Read points.
  */
  rewind( fp );
  for( i = 0; i < numPoint; ++i )
  {
    nitem = fread( pbuf, POINT_SIZE, 1, fp );
    if( nitem != 1 )
    {
      fprintf( stderr,
               "readESD: unable to read data point %d in file %s\n",
               i + 1,
               filename );
      fclose( fp );
      chromatData->numPoint = 0;
      chromatData->numBase  = 0;
      *status = 1;
      return( chromatData );
    }

    cptr = pbuf + 4;
    chromatData->trace[baseChar[(int)baseMap[0]]][i] = (FLOAT)inSwpFloatL( cptr );

    cptr += 4;
    chromatData->trace[baseChar[(int)baseMap[1]]][i] = (FLOAT)inSwpFloatL( cptr );

    cptr += 4;
    chromatData->trace[baseChar[(int)baseMap[2]]][i] = (FLOAT)inSwpFloatL( cptr );

    cptr += 4;
    chromatData->trace[baseChar[(int)baseMap[3]]][i] = (FLOAT)inSwpFloatL( cptr );
  }

  /*
  ** Find minimum and maximum trace values.
  */
  findTraceExtrema( chromatData );

  /*
  ** Close file.
  */
  fclose( fp );

  if( option->verboseOption && option->verboseLevel >= 16 )
  {
    fprintf( stderr,
             "readESD: machine name:     %s\n",
             chromatData->machineName );
    fprintf( stderr,
             "readESD: gel name:         %s\n",
             chromatData->gelName );
    fprintf( stderr,
             "readESD: sample name:      %s\n",
             chromatData->sampleName );
    fprintf( stderr,
             "readESD: primer ID:        %s\n",
             chromatData->primerID );
    fprintf( stderr,
             "readESD: well ID:          %s\n",
              extHeader.wellID );
    fprintf( stderr,
             "readESD: comment:          %s\n",
             chromatData->comment );
    fprintf( stderr,
             "readESD: number of scans:  %d\n",
             chromatData->numPoint );
    fprintf( stderr,
             "readESD: base order:       %c %c %c %c\n",
             baseMap[0],
             baseMap[1],
             baseMap[2],
             baseMap[3] );
  }

  *status = 0;

  return( chromatData );
}

