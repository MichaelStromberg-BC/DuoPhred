/** chromat2phred.c **/

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
#ifndef ANSI_C
#include <memory.h>
#endif
#include "phred.h"

/*
** Move required elements of chromatData to
** phredData.
*/
#ifdef ANSI_C
PhredData *chromat2phred( ChromatData *chromatData, Option *option )
#else
PhredData *chromat2phred( chromatData, option )
ChromatData *chromatData;
Option *option;
#endif
{
  int i, j;
  int flag;
  int numByte;
  int parMatchFlag;
  char message[2*PHRED_PATH_MAX+2*PHRED_MAX_STRING_LEN];
  PhredData *phredData;

  /*
  ** Create phredData structure.
  */
  numByte = sizeof( PhredData );
  phredData = (PhredData *)ourMalloc( numByte );

  /*
  ** Initialize structure.
  */
  phredData->numPoint = 0;
  for( j = 0; j < 4; ++j )
  {
    phredData->trace[j] = NULL;
  }
  for( j = 0; j < 2; ++j )
  {
    phredData->numBase[j] = 0;
    phredData->base[j] = NULL;
	phredData->baseUncalled[j] = NULL;
    phredData->baseLoc[j] = NULL;
	phredData->baseLocUncalled[j] = NULL;
    phredData->baseQual[j] = NULL;
	phredData->baseQualUncalled[j] = NULL;
  }
  phredData->numDsc = 0;
  phredData->baseDsc = NULL;
  phredData->qualIndex = NULL;
  phredData->polyData.called_nuc             = NULL;
  phredData->polyData.called_loc             = NULL;
  phredData->polyData.called_area            = NULL;
  phredData->polyData.called_relative_area   = NULL;
  phredData->polyData.uncalled_nuc           = NULL;
  phredData->polyData.uncalled_loc           = NULL;
  phredData->polyData.uncalled_area          = NULL;
  phredData->polyData.uncalled_relative_area = NULL;
  phredData->polyData.A_tr_val               = NULL;
  phredData->polyData.C_tr_val               = NULL;
  phredData->polyData.G_tr_val               = NULL;
  phredData->polyData.T_tr_val               = NULL;

  phredData->lftPhredTrim                    = -1;
  phredData->rhtPhredTrim                    = -1;
  phredData->trimSetValue                    = 0.0;

  /*
  ** Copy trace data to phredData.
  */
  phredData->numPoint = chromatData->numPoint;
  if( chromatData->numPoint > 0 )
  {
    numByte = chromatData->numPoint * sizeof( FLOAT );
    for( j = 0; j < 4; ++j )
    {
      phredData->trace[j] = (FLOAT *)ourMalloc( numByte );
      memcpy( phredData->trace[j], chromatData->trace[j], numByte );
    }
  }

  /*
  ** Copy bases to phredData.
  */
  phredData->numBase[FIL] = chromatData->numBase;
  if( chromatData->numBase > 0 )
  {
    numByte = chromatData->numBase * sizeof( char );
    phredData->base[FIL] = (char *)ourMalloc( numByte );
    memcpy( phredData->base[FIL], chromatData->base, numByte );
  }

  /*
  ** Copy base locations to phredData.
  */
  if( chromatData->numBase > 0 )
  {
    numByte = chromatData->numBase * sizeof( int );
    phredData->baseLoc[FIL] = (int *)ourMalloc( numByte );
    memcpy( phredData->baseLoc[FIL], chromatData->baseLoc, numByte );
  }

  /*
  ** Copy base quality to phredData.
  */
  if( chromatData->numBase > 0 )
  {
    numByte = chromatData->numBase * sizeof( int );
    phredData->baseQual[FIL] = (int *)ourMalloc( numByte );
    memcpy( phredData->baseQual[FIL], chromatData->baseQual, numByte );
  }

  /*
  ** Allocate memory.
  */
  numByte = chromatData->numBase * sizeof( FLOAT );
  phredData->qualIndex = (FLOAT *)ourMalloc( numByte );

  /*
  ** Copy information.
  */
  phredData->primerLoc     = chromatData->primerLoc;
  phredData->avgSpacing    = chromatData->avgSpacing;
  phredData->minTraceValue = chromatData->minTraceValue;
  phredData->maxTraceValue = chromatData->maxTraceValue;
  phredData->fileType      = chromatData->fileType;
  phredData->laneNumber    = chromatData->laneNumber;

  for( j = 0; j < 4; ++j )
  {
    phredData->signalStrength[j] = chromatData->signalStrength[j];
  }
  pstrcpy( phredData->fileName,    chromatData->fileName,    PHRED_PATH_MAX );
  pstrcpy( phredData->comment,     chromatData->comment,     PHRED_MAX_STRING_LEN );
  pstrcpy( phredData->convProg,    chromatData->convProg,    PHRED_MAX_STRING_LEN );
  pstrcpy( phredData->gelName,     chromatData->gelName,     PHRED_MAX_STRING_LEN );
  pstrcpy( phredData->machineName, chromatData->machineName, PHRED_MAX_STRING_LEN );
  pstrcpy( phredData->primerID,    chromatData->primerID,    PHRED_MAX_STRING_LEN );
  pstrcpy( phredData->processing,  chromatData->processing,  PHRED_MAX_STRING_LEN );
  pstrcpy( phredData->reTracker,   chromatData->reTracker,   PHRED_MAX_STRING_LEN );
  pstrcpy( phredData->sampleName,  chromatData->sampleName,  PHRED_MAX_STRING_LEN );
  memcpy( phredData->thumbPrint,  chromatData->thumbPrint, 10 );

  /*
  ** Initialize some values.
  */
  phredData->leftTrimPoint     = 0;
  phredData->rghtTrimPoint     = 0;
  phredData->compressSplitFlag = option->compressSplitFlag;
  phredData->primerQVFlag      = option->primerQVFlag;

  /*
  ** Try to work around bogus primer id string:
  ** I see cases where the id string consists of
  ** only non-printing characters. Interpret this
  ** as a NULL string.
  */
  flag = 0;
  for( i = 0; i < strlen( phredData->primerID ); ++i )
  {
    if( phredData->primerID[i] > 32 )
    {
      flag = 1;
      break;
    }
  }
  if( flag == 0 )
  {
    phredData->primerID[0] = '\0';
  }

  /*
  ** Determine chemistry, if requested.
  */
  parMatchFlag           = 0;
  phredData->chemType    = UNKNOWN_CHEM;
  phredData->dyeType     = UNKNOWN_DYE;
  phredData->machineType = UNKNOWN_MACHINE;
  for( j = 0; j < option->parFileData.numChem; ++j )
  {
    if( strcmp( phredData->primerID, option->parFileData.chemList[j].name ) == 0 )
    {
      phredData->chemType    = option->parFileData.chemList[j].chemType;
      phredData->dyeType     = option->parFileData.chemList[j].dyeType;
      phredData->machineType = option->parFileData.chemList[j].machineType;
      parMatchFlag           = 1;
      break;
    }
  }

  /*
  ** Set machine type.
  */
  if( phredData->machineType == ABI_373_377 )
  {
    pstrcpy( phredData->source, "ABI 373A or 377", PHRED_MAX_STRING_LEN );
  }
  else
  if( phredData->machineType == MOLDYN_MEGABACE )
  {
    pstrcpy( phredData->source, "MegaBACE", PHRED_MAX_STRING_LEN );
  }
  else
  if( phredData->machineType == ABI_3700 )
  {
    pstrcpy( phredData->source, "ABI 3700", PHRED_MAX_STRING_LEN );
  }
  else
  if( phredData->machineType == LI_COR_4000 )
  {
    pstrcpy( phredData->source, "LI-COR 4000", PHRED_MAX_STRING_LEN );
  }
  else
  if( phredData->machineType == UNKNOWN_MACHINE )
  {
    pstrcpy( phredData->source, "unknown", PHRED_MAX_STRING_LEN );
  }

  /*
  ** Diagnostics.
  */
  if( option->verboseOption && option->verboseLevel >= 1 )
  {
    fprintf( stderr,
             "chromat2phred: chemistry type:      " );
    if( phredData->chemType == PRIMER_CHEM )
    {
      fprintf( stderr,
               "dye primer\n" );
    }
    else
    if( phredData->chemType == TERMINATOR_CHEM )
    {
      fprintf( stderr,
               "dye terminator\n" );
    }
    else
    if( phredData->chemType == UNKNOWN_CHEM )
    {
      fprintf( stderr,
               "unknown chemistry\n" );
    }

    fprintf( stderr,
             "chromat2phred: dye type:            " );
    if( phredData->dyeType == RHODAMINE_DYE )
    {
      fprintf( stderr,
               "rhodamine\n" );
    }
    else
    if( phredData->dyeType == BIG_DYE_DYE )
    {
      fprintf( stderr,
               "big dye\n" );
    }
    else
    if( phredData->dyeType == D_RHODAMINE_DYE )
    {
      fprintf( stderr,
               "d-rhodamine\n" );
    }
    else
    if( phredData->dyeType == ENERGY_TRANSFER_DYE )
    {
      fprintf( stderr,
               "energy transfer\n" );
    }
    else
    if( phredData->dyeType == BODIPY_DYE )
    {
      fprintf( stderr,
               "bodipy\n" );
    }
    else
    if( phredData->dyeType == UNKNOWN_DYE )
    {
      fprintf( stderr,
               "unknown dye\n" );
    }

    fprintf( stderr,
             "chromat2phred: sequencing machine:  " );
    if( phredData->machineType == ABI_373_377 )
    {
      fprintf( stderr,
               "ABI 373 or ABI 377\n" );
    }
    else
    if( phredData->machineType == MOLDYN_MEGABACE )
    {
      fprintf( stderr,
               "MegaBACE\n" );
    }
    else
    if( phredData->machineType == ABI_3700 )
    {
      fprintf( stderr,
               "ABI 3700\n" );
    }
    else
    if( phredData->machineType == LI_COR_4000 )
    {
      fprintf( stderr,
               "LI-COR\n" );
    }
    else
    if( phredData->machineType == UNKNOWN_MACHINE )
    {
      fprintf( stderr,
               "unknown machine\n" );
    }

  }

  if( option->verboseOption && option->verboseLevel >= 32 )
  {
    fprintf( stderr,
             "chromat2phred: minTraceValue:       %f\n",
             phredData->minTraceValue );
    fprintf( stderr,
             "chromat2phred: maxTraceValue:       %f\n",
             phredData->maxTraceValue );
  }

  /*
  ** Emit an appropriate error message if
  ** we could not match the primer ID string.
  */
  if( parMatchFlag == 0 )
  {
    if( option->parFileName[0] == '\0' )
    {
      /*
      ** PHRED_PARAMETER_FILE environment variable unset.
      */
      sprintf( message,
               "warning: 'PHRED_PARAMETER_FILE' environment variable not set:\n" );
      fprintf( stderr, "%s", message );
      if( option->log )
      {
        writeLog( message );
      }

      sprintf( message,
             "          unable to identify chemistry and dye\n" );
      fprintf( stderr, "%s", message );
      if( option->log )
      {
        writeLog( message );
      }

      sprintf( message,
               "          type `phred -doc' for more information\n" );
      fprintf( stderr, "%s", message );
      if( option->log )
      {
        writeLog( message );
      }
    }
    else
    if( option->parFileData.successfulRead == 0 )
    {
      /*
      ** Unable to read phred parameter file. Either it
      ** doesn't exist or it is unreadable/parsable.
      */
      sprintf( message,
               "warning: unable to read parameter file: %.*s\n",
               PHRED_PATH_MAX,
               option->parFileName );
      fprintf( stderr, "%s", message );
      if( option->log )
      {
        writeLog( message );
      }

      sprintf( message,
               "         unable to identify chemistry and dye\n" );
      fprintf( stderr, "%s", message );
      if( option->log )
      {
        writeLog( message );
      }

      sprintf( message,
               "         type `phred -doc' for more information\n" );
      fprintf( stderr, "%s", message );
      if( option->log )
      {
        writeLog( message );
      }
    }
    else
    if( phredData->numPoint > 0 &&
        option->parFileData.successfulRead == 1 )
    {
      /*
      ** Unable to find a primer ID string match in the
      ** phred parameter file entries.
      */
      if( phredData->primerID[0] != '\0' )
      {
        sprintf( message,
                 "    unknown chemistry (%.*s) in chromat %.*s\n",
                 PHRED_MAX_STRING_LEN,
                 phredData->primerID,
                 PHRED_PATH_MAX,
                 option->inFileName[option->curInFile] );
        fprintf( stderr, "%s", message );
        if( option->log )
        {
          writeLog( message );
        }

        sprintf( message,
                 "    add a line of the form\n" );
        fprintf( stderr, "%s", message );
        if( option->log )
        {
          writeLog( message );
        }

        sprintf( message,
                 "    \"%.*s\"    <chemistry>      <dye type>      <machine type>\n",
                 PHRED_MAX_STRING_LEN,
                 phredData->primerID );
        fprintf( stderr, "%s", message );
        if( option->log )
        {
          writeLog( message );
        }

        sprintf( message,
                 "    to the file %.*s\n",
                 PHRED_PATH_MAX,
                 option->parFileName );
        fprintf( stderr, "%s", message );
        if( option->log )
        {
          writeLog( message );
        }

        sprintf( message,
                 "    type `phred -doc' for more information\n" );
        fprintf( stderr, "%s", message );
        if( option->log )
        {
          writeLog( message );
        }
      }
      else
      {
        sprintf( message,
                 "    no dye primer ID in chromat %.*s\n",
                 PHRED_PATH_MAX,
                 option->inFileName[option->curInFile] );
        fprintf( stderr, "%s", message );
        if( option->log )
        {
          writeLog( message );
        }
      }
    }
  }

  return( phredData );
}

