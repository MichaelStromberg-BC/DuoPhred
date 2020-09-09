/** readParam.c **/

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
#include "phred.h"

/*
**
** See helpParam.c for parameters and descriptions.
**
** Input files are specified on the command line.
**
** Add command line parameters where noted below.
**
** Note: initialize options in initParam.c
**
*/

#ifdef ANSI_C
static char *getArg( char *arg )
#else
static char *getArg( arg )
char *arg;
#endif
{
  if( arg[0] == '-' )
  {
    return( NULL );
  }
  return( arg );
}



#ifdef ANSI_C
int readParam( int argc, char **argv, Option *option )
#else
int readParam( argc, argv, option )
int argc;
char **argv;
Option *option;
#endif
{
  int iarg;
  char *string;

  /*
  ** read through command line parameters
  */
  iarg = 1;
  while( iarg < argc )
  {
    /*
    ** options begin with a dash
    */
    if( strcmp( argv[iarg], "-view" ) == 0 )
    {
      option->edit = 1;
    }
    else
    if( strcmp( argv[iarg], "-s" ) == 0 )
    {
      if( iarg+1 >= argc )
      {
        option->writeSeq = 1;
        option->sOption = 1;
        option->writeSeqName = NULL;
      }
      else
      {
        string = getArg( argv[iarg+1] );
        option->writeSeq = 1;
        option->writeSeqName = string;
        if( string != NULL )
        {
          ++iarg;
          option->sOption = 2;
        }
        else
        {
          option->sOption = 1;
        }
      }
    }
    else
    if( strcmp( argv[iarg], "-sa" ) == 0 )
    {
      option->writeSeq = 1;
      option->saOption = 1;
      ++iarg;
      if( iarg >= argc )
      {
        fprintf( stderr, "-sa: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      if( ( string = getArg( argv[iarg] ) ) == NULL )
      {
        fprintf( stderr, "-sa: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      option->writeSeqName = string;
    }
    else
    if( strcmp( argv[iarg], "-sd" ) == 0 )
    {
      option->writeSeq = 1;
      option->sdOption = 1;
      ++iarg;
      if( iarg >= argc )
      {
        fprintf( stderr, "-sd: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      if( ( string = getArg( argv[iarg] ) ) == NULL )
      {
        fprintf( stderr, "-sd: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      option->seqDirName = string;
    }
    else
    if( strcmp( argv[iarg], "-q" ) == 0 )
    {
      if( iarg+1 >= argc )
      {
        option->writeQual = 1;
        option->qOption = 1;
        option->writeQualName = NULL;
      }
      else
      {
        string = getArg( argv[iarg+1] );
        option->writeQual = 1;
        option->writeQualName = string;
        if( string != NULL )
        {
          option->qOption = 2;
          ++iarg;
        }
        else
        {
          option->qOption = 1;
        }
      }
    }
    else
    if( strcmp( argv[iarg], "-qa" ) == 0 )
    {
      option->writeQual = 1;
      option->qaOption = 1;
      ++iarg;
      if( iarg >= argc )
      {
        fprintf( stderr, "-qa: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      if( ( string = getArg( argv[iarg] ) ) == NULL )
      {
        fprintf( stderr, "-qa: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      option->writeQualName = string;
    }
    else
    if( strcmp( argv[iarg], "-qd" ) == 0 )
    {
      option->writeQual = 1;
      option->qdOption = 1;
      ++iarg;
      if( iarg >= argc )
      {
        fprintf( stderr, "-qd: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      if( ( string = getArg( argv[iarg] ) ) == NULL )
      {
        fprintf( stderr, "-qd: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      option->qualDirName = string;
    }
    else
    if( strcmp( argv[iarg], "-cp" ) == 0 )
    {
      option->scfPrecOption = 1;
      ++iarg;
      if( iarg >= argc )
      {
        fprintf( stderr, "-cp: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      if( ( string = getArg( argv[iarg] ) ) == NULL )
      {
        fprintf( stderr, "-cp: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      option->scfPrecision = atoi( string );
    }
    else
    if( strcmp( argv[iarg], "-cv" ) == 0 )
    {
      ++iarg;
      if( iarg >= argc )
      {
        fprintf( stderr, "-cv: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      if( ( string = getArg( argv[iarg] ) ) == NULL )
      {
        fprintf( stderr, "-cv: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      option->scfVersion = atoi( string );
    }
    else
    if( strcmp( argv[iarg], "-cs" ) == 0 )
    {
      option->scaleSCFTraceOption = 1;
    }
    else
    if( strcmp( argv[iarg], "-c" ) == 0 )
    {
      if( iarg+1 >= argc )
      {
        option->writeScf = 1;
        option->writeScfName = NULL;
      }
      else
      {
        string = getArg( argv[iarg+1] );
        option->writeScf = 1;
        option->writeScfName = string;
        if( string != NULL )
        {
          ++iarg;
        }
      }
    }
    else
    if( strcmp( argv[iarg], "-cd" ) == 0 )
    {
      option->writeScf = 1;
      option->scfDir = 1;
      ++iarg;
      if( iarg >= argc )
      {
        fprintf( stderr, "-cd: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      if( ( string = getArg( argv[iarg] ) ) == NULL )
      {
        fprintf( stderr, "-cd: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      option->scfDirName = string;
    }
    else
    if( strcmp( argv[iarg], "-p" ) == 0 )
    {
      option->writePhd = 1;
      option->pOption = 1;
      if( iarg+1 >= argc )
      {
        option->writePhdName = NULL;
      }
      else
      {
        string = getArg( argv[iarg+1] );
        option->writePhdName = string;
        if( string != NULL )
        {
          ++iarg;
        }
      }
    }
    else
    if( strcmp( argv[iarg], "-pd" ) == 0 )
    {
      option->writePhd = 1;
      option->pdOption = 1;
      ++iarg;
      if( iarg >= argc )
      {
        fprintf( stderr, "-pd: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      if( ( string = getArg( argv[iarg] ) ) == NULL )
      {
        fprintf( stderr, "-pd: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      option->phdDirName = string;
    }
    else
	if( strcmp( argv[iarg], "-x" ) == 0 )
    {
		option->writeXmlData =1;
      option->xOption = 1;
      if( iarg+1 >= argc )
      {
		  option->writeXmlDataName = NULL;
      }
      else
      {
        string = getArg( argv[iarg+1] );
        option->writeXmlDataName = string;
        if( string != NULL )
        {
          ++iarg;
        }
      }
    }
    else
	if( strcmp( argv[iarg], "-xd" ) == 0 )
    {
	  option->writeXmlData = 1;
      option->xdOption = 1;
      ++iarg;
      if( iarg >= argc )
      {
        fprintf( stderr, "-xd: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      if( ( string = getArg( argv[iarg] ) ) == NULL )
      {
        fprintf( stderr, "-xd: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
	  option->xmlDataDirName = string;
	  //printf("MPS: Finished the XML parameter allocation.\n");
    }
    else
    if( strcmp( argv[iarg], "-bottom" ) == 0 )
    {
      option->bottom = 1;
    }
    else
    if( strcmp( argv[iarg], "-trim" ) == 0 )
    {
      if( option->trim != 0 )
      {
        fprintf( stderr,
                 "-trim: override %s\n",
                 option->trim == 1 ? "-trim" : "-trim_alt" );
      }
      option->trim = 1;
      ++iarg;
      if( iarg >= argc )
      {
        fprintf( stderr, "-trim: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      if( ( string = getArg( argv[iarg] ) ) == NULL )
      {
        fprintf( stderr, "-trim: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      option->enzName = string;
    }
    else
    if( strcmp( argv[iarg], "-trim_alt" ) == 0 )
    {
      if( option->trim != 0 )
      {
        fprintf( stderr,
                 "-trim_alt: override %s\n",
                 option->trim == 1 ? "-trim" : "-trim_alt" );
      }
      option->trim = 2;
      ++iarg;
      if( iarg >= argc )
      {
        fprintf( stderr, "-trim_alt: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      if( ( string = getArg( argv[iarg] ) ) == NULL )
      {
        fprintf( stderr, "-trim_alt: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      option->enzName = string;
    }
    else
    if( strcmp( argv[iarg], "-nocall" ) == 0 )
    {
      option->call = 0;
    }
    else
    if( strcmp( argv[iarg], "-st" ) == 0 )
    {
      ++iarg;
      if( iarg >= argc )
      {
        fprintf( stderr, "-st: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      if( ( string = getArg( argv[iarg] ) ) == NULL )
      {
        fprintf( stderr, "-st: argument missing\n" );
        option->errorFlag = 1;
        continue;     
      }
      if( strcmp( string, "fasta" ) == 0 )
      {
        option->seqType = 0;
      }
      else
      if( strcmp( string, "xbap" ) == 0 )
      {
        option->seqType = 1;
      }
      else
      {
        fprintf( stderr, "-st: bad argument\n" );
        option->errorFlag = 1;
        continue;
      }
    }
    else
    if( strcmp( argv[iarg], "-qt" ) == 0 )
    {
      ++iarg;
      if( iarg >= argc )
      {
        fprintf( stderr, "-qt: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      if( ( string = getArg( argv[iarg] ) ) == NULL )
      {
        fprintf( stderr, "-qt: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      if( strcmp( string, "fasta" ) == 0 )
      {
        option->qualType = 0;
      }
      else
      if( strcmp( string, "xbap" ) == 0 )
      {
        option->qualType = 1;
      }
      else
      if( strcmp( string, "mix" ) == 0 )
      {
        option->qualType = 2;
      }
      else
      {
        fprintf( stderr, "-qt: bad argument\n" );
        option->errorFlag = 1;
        continue;
      }
    }
    else
    if( strcmp( argv[iarg], "-help" ) == 0 ||
        strcmp( argv[iarg], "-h" ) == 0 )
    {
      helpParam();
    }
    else
    if( strcmp( argv[iarg], "-mag" ) == 0 )
    {
      option->magnification = 30;
      ++iarg;
      if( iarg >= argc )
      {
        fprintf( stderr, "-mag: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      if( ( string = getArg( argv[iarg] ) ) == NULL )
      {
        fprintf( stderr, "-mag: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      option->magnification = atoi( string );
    }
    else
    if( strcmp( argv[iarg], "-raw" ) == 0 )
    {
      option->seq = 1;
      ++iarg;
      if( iarg >= argc )
      {
        fprintf( stderr, "-raw: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      if( ( string = getArg( argv[iarg] ) ) == NULL )
      {
        fprintf( stderr, "-raw: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      option->seqName = string;
    }
    else
    if( strcmp( argv[iarg], "-bn" ) == 0 )
    {
      option->baseNumber = 1;
      ++iarg;
      if( iarg >= argc )
      {
        fprintf( stderr, "-bn: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      if( ( string = getArg( argv[iarg] ) ) == NULL )
      {
        fprintf( stderr, "-bn: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      option->baseNumber = atoi( string );
    }
    else
    if( strcmp( argv[iarg], "-loc" ) == 0 )
    {
      option->subSeq = 1;
      ++iarg;
      if( iarg >= argc )
      {
        fprintf( stderr, "-loc: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      if( ( string = getArg( argv[iarg] ) ) == NULL )
      {
        fprintf( stderr, "-loc: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      option->subSequence = string;
    }
    else
    if( strcmp( argv[iarg], "-log" ) == 0 )
    {
      option->log = 1;
    }
    else
    if( strcmp( argv[iarg], "-if" ) == 0 )
    {
      option->ifOption = 1;
      ++iarg;
      if( iarg >= argc )
      {
        fprintf( stderr, "-if: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      if( ( string = getArg( argv[iarg] ) ) == NULL )
      {
        fprintf( stderr, "-if: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      if( option->idOption || readFOF( argc, string, option ) == ERROR )
      {
        option->errorFlag = 1;
        continue;
      }
    }
    else
    if( strcmp( argv[iarg], "-id" ) == 0 )
    {
      option->idOption = 1;
      ++iarg;
      if( iarg >= argc )
      {
        fprintf( stderr, "-id: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      if( ( string = getArg( argv[iarg] ) ) == NULL )
      {
        fprintf( stderr, "-id: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      if( option->ifOption || readDIR( argc, string, option ) == ERROR )
      {
        option->errorFlag = 1;
        continue;
      }
    }
    else
    if( strcmp( argv[iarg], "-sync" ) == 0 )
    {
      /*
      ** Pass through to X code by acknowledging this option
      ** but taking no action.
      */
    }
    else
    if( strcmp( argv[iarg], "-V" ) == 0 )
    {
      printf( "\n" );
      printf( "  phred version: %s\n", getVersion() );
      printf( "\n" );
    }
    else
    if( strcmp( argv[iarg], "-qr" ) == 0 )
    {
      ++iarg;
      if( iarg >= argc )
      {
        fprintf( stderr, "-qr: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      if( ( string = getArg( argv[iarg] ) ) == NULL )
      {
        fprintf( stderr, "-qr: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      option->qualReportFileName = string;
      option->qualReport = 1;
    }
    else
    if( strcmp( argv[iarg], "-d" ) == 0 )
    {
      option->writePolyData = 1;
      option->dOption = 1;
      if( iarg+1 >= argc )
      {
        option->writePolyDataName = NULL;
      }
      else
      {
        string = getArg( argv[iarg+1] );
        option->writePolyDataName = string;
        if( string != NULL )
        {
          ++iarg;
        }
      }
    }
    else
    if( strcmp( argv[iarg], "-dd" ) == 0 )
    {
      option->writePolyData = 1;
      option->ddOption = 1;
      ++iarg;
      if( iarg >= argc )
      {
        fprintf( stderr, "-dd: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      if( ( string = getArg( argv[iarg] ) ) == NULL )
      {
        fprintf( stderr, "-dd: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      option->polyDataDirName = string;
    }
    else
    if( strcmp( argv[iarg], "-xm" ) == 0 )
    {
      option->xmOption = 1;
      ++iarg;
      if( iarg >= argc )
      {
        fprintf( stderr, "-xm: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      if( ( string = getArg( argv[iarg] ) ) == NULL )
      {
        fprintf( stderr, "-xm: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      option->xmFileName = string;
    }
    else
    if( strcmp( argv[iarg], "-zd" ) == 0 )
    {
      option->filCompressExeDirOption = 1;
      ++iarg;
      if( iarg >= argc )
      {
        fprintf( stderr, "-zd: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      if( ( string = getArg( argv[iarg] ) ) == NULL )
      {
        fprintf( stderr, "-zd: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      option->filCompressExeDir = string;
    }
    else
    if( strcmp( argv[iarg], "-zt" ) == 0 )
    {
      option->filCompressTmpDirOption = 1;
      ++iarg;
      if( iarg >= argc )
      {
        fprintf( stderr, "-zt: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      if( ( string = getArg( argv[iarg] ) ) == NULL )
      {
        fprintf( stderr, "-zt: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      option->filCompressTmpDir = string;
    }
    else
    if( strcmp( argv[iarg], "-nonorm" ) == 0 )
    {
      option->normalize = 0;
    }
    else
    if( strcmp( argv[iarg], "-nosplit" ) == 0 )
    {
      option->compressSplitFlag = 0;
    }
    else
    if( strcmp( argv[iarg], "-nocmpqv" ) == 0 )
    {
      option->primerQVFlag = 0;
    }
    else
    if( strcmp( argv[iarg], "-ceilqv" ) == 0 )
    {
      option->qualityValueCeilingOption = 1;
      ++iarg;
      if( iarg >= argc )
      {
        fprintf( stderr, "-ceilqv: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      if( ( string = getArg( argv[iarg] ) ) == NULL )
      {
        fprintf( stderr, "-ceilqv: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      option->qualityValueCeiling = atoi( string );
    }
    else
    if( strcmp( argv[iarg], "-beg_pred" ) == 0 )
    {
      option->beginPeakPredictionOption = 1;
      ++iarg;
      if( iarg >= argc )
      {
        fprintf( stderr, "-beg_pred: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      if( ( string = getArg( argv[iarg] ) ) == NULL )
      {
        fprintf( stderr, "-beg_pred: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      option->beginPeakPredictionPoint = atoi( string );
    }
    else
    if( strcmp( argv[iarg], "-tags" ) == 0 )
    {
      option->tagOption = 1;
    }
    else
    if( strcmp( argv[iarg], "-doc" ) == 0 )
    {
      showDoc();
    }
    else
    if( strcmp( argv[iarg], "-v" ) == 0 )
    {
      option->verboseOption = 1;
      ++iarg;
      if( iarg >= argc )
      {
        fprintf( stderr, "-v: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      if( ( string = getArg( argv[iarg] ) ) == NULL )
      {
        fprintf( stderr, "-v: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      option->verboseLevel = atoi( string );
    }
    else
    if( strcmp( argv[iarg], "-trim_cutoff" ) == 0 )
    {
      option->trimSetOption = 1;
      ++iarg;
      if( iarg >= argc )
      {
        fprintf( stderr, "-trim_cutoff: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      if( ( string = getArg( argv[iarg] ) ) == NULL )
      {
        fprintf( stderr, "-trim_cutoff: argument missing\n" );
        option->errorFlag = 1;
        continue;
      }
      option->trimSetValue = atof( string );
    }
    else
    if( strcmp( argv[iarg], "-trim_fasta" ) == 0 )
    {
      option->trimFastaData = 1;
    }
    else
    if( strcmp( argv[iarg], "-trim_scf" ) == 0 )
    {
      option->trimSCFData = 1;
    }
    else
    if( strcmp( argv[iarg], "-trim_phd" ) == 0 )
    {
      option->trimPHDData = 1;
    }
    else
    if( strcmp( argv[iarg], "-trim_out" ) == 0 )
    {
      option->trimFastaData = 1;
      option->trimSCFData   = 1;
      option->trimPHDData   = 1;
    }
    
/*
** add options above here
*/
    else
    if( argv[iarg][0] == '-' )
    {
      fprintf( stderr, "unknown option %s\n", argv[iarg] );
      option->errorFlag = 1;
    }
    else
    {
      /*
      ** input filenames
      */
      option->inFileName[option->numInFile] = argv[iarg];
      ++option->numInFile;
    }
    ++iarg;
  }

  return( OK );
}
