/** initParam.c **/

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

#ifdef ANSI_C
int initParam( Option *option, int argc )
#else
int initParam( option, argc )
Option *option;
int argc;
#endif
{

  /*
  ** initialize program parameters to defaults
  */
  option->edit                       = 0;
  option->editNumber                 = -1;
  option->writeSeq                   = 0;
  option->sOption                    = 0;
  option->saOption                   = 0;
  option->sdOption                   = 0;
  option->writeSeqName               = NULL;
  option->seqDirName                 = NULL;
  option->writeQual                  = 0;
  option->qOption                    = 0;
  option->qaOption                   = 0;
  option->qdOption                   = 0;
  option->writeQualName              = NULL;
  option->qualDirName                = NULL;
  option->writeScf                   = 0;
  option->writeScfName               = NULL;
  option->writePolyData              = 0;
  option->writeXmlData				 = 0;
  option->dOption                    = 0;
  option->ddOption                   = 0;
  option->xOption                    = 0;
  option->xdOption                   = 0;
  option->writePolyDataName          = NULL;
  option->writeXmlDataName           = NULL;
  option->polyDataDirName            = NULL;
  option->xmlDataDirName            = NULL;
  option->scfDir                     = 0;
  option->scfDirName                 = NULL;
  option->writePhd                   = 0;
  option->pOption                    = 0;
  option->pdOption                   = 0;
  option->writePhdName               = NULL;
  option->phdDirName                 = NULL;
  option->bottom                     = 0;
  option->trim                       = 0;
  option->enzName                    = NULL;
  option->call                       = 1;
  option->normalize                  = 1;
  option->seqType                    = 0;
  option->qualType                   = 0;
  option->inType                     = 0;
  option->ifOption                   = 0;
  option->idOption                   = 0;
  option->numInFile                  = 0;
  option->diag                       = 0;
  option->log                        = 0;
  option->errorFlag                  = 0;
  option->seq                        = 0;
  option->magnification              = 0;
  option->subSeq                     = 0;
  option->subSequence                = NULL;
  option->baseNumber                 = 0;
  option->seqName                    = NULL;
  option->newpred                    = 1;
  option->qualReport                 = 0;
  option->qualReportFileName         = NULL;
  option->scfPrecOption              = 0;
  option->scfPrecision               = 2;
  option->scfVersion                 = 2;
  option->scaleSCFTraceOption        = 0;
  option->xmOption                   = 0;
  option->xmFileName                 = NULL;
  option->filCompressExeDirOption    = 0;
  option->filCompressExeDir          = NULL;
  option->filCompressTmpDirOption    = 0;
  option->filCompressTmpDir          = NULL;
  option->parFileName[0]             = '\0';
  option->parFileData.numChem        = 0;
  option->parFileData.chemList       = NULL;
  option->parFileData.successfulRead = 0;
  option->compressSplitFlag          = 1;
  option->primerQVFlag               = 1;
  option->qualityValueCeilingOption  = 0;
  option->qualityValueCeiling        = 100;
  option->beginPeakPredictionOption  = 0;
  option->beginPeakPredictionPoint   = 0;
  option->tagOption                  = 0;
  option->verboseOption              = 0;
  option->verboseLevel               = 0;
  option->trimSetOption              = 0;
  option->trimSetValue               = 0;
  option->trimFastaData              = 0;
  option->trimSCFData                = 0;
  option->trimPHDData                = 0;


  if( ( option->inFileName = (char **)malloc( argc * sizeof( char * ) ) ) == NULL )
  {
    fprintf( stderr, "initParam: unable to allocate memory\n" );
    return( -1 );
  }

  initFOF();
  initDIR();

  return( 0 );

}
