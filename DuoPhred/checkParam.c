/** checkParam.c **/

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
** Check option consistency.
*/

#ifdef ANSI_C
int checkParam( Option *option )
#else
int checkParam( option )
Option *option;
#endif
{
  /*
  ** consistency tests
  */
  if( option->ifOption && option->idOption )
  {
    fprintf( stderr,
             "use one of -if and -id\n" );
    if( option->log )
    {
      writeLog( "use one of -if and -id\n" );
    }
  }
  if( ( option->sOption  && option->saOption ) ||
      ( option->sOption  && option->sdOption ) ||
      ( option->saOption && option->sdOption ) )
  {
    fprintf( stderr,
             "use one of -s, -sa, and -sd\n" );
    if( option->log )
    {
      writeLog( "use one of -s, -sa, and -sd\n" );
    }
    option->errorFlag = 1;
  }
  if( ( option->qOption  && option->qaOption ) ||
      ( option->qOption  && option->qdOption ) ||
      ( option->qaOption && option->qdOption ) )
  {
    fprintf( stderr,
             "use one of -q, -qa, and -qd\n" );
    if( option->log )
    {
      writeLog( "use one of -q, -qa, and -qd\n" );
    }
    option->errorFlag = 1;
  }

  if( option->numInFile > 1 && option->sOption == 2 )
  {
    fprintf( stderr,
             "use default sequence filename with more than one input file\n" );
    if( option->log )
    {
      writeLog( "use default sequence filename with more than one input file\n" );
    }
    option->errorFlag = 1;
  }
  if( option->numInFile > 1 && option->qOption == 2 )
  {
    fprintf( stderr,
             "use default quality filename with more than one input file\n" );
    if( option->log )
    {
      writeLog( "use default quality filename with more than one input file\n" );
    }
    option->errorFlag = 1;
  }
  if( option->seqType == 1 && option->saOption )
  {
    fprintf( stderr,
             "-st xbap excludes -sa\n" );
    if( option->log )
    {
      writeLog( "-st xbap excludes -sa\n" );
    }
    option->errorFlag = 1;
  }
  if( option->numInFile > 1 && option->writeScfName != NULL )
  {
    fprintf( stderr,
             "use -cd with more than one input file\n" );
    if( option->log )
    {
      writeLog( "use -cd with more than one input file\n" );
    }
    option->errorFlag = 1;
  }
  if( option->numInFile == 0 && option->edit == 0 )
  {
    fprintf( stderr,
             "no input files specified\n" );
    if( option->log )
    {
      writeLog( "no input files specified\n" );
    }
    option->errorFlag = 1;
  }
  if( option->numInFile > 1 && option->seq )
  {
    fprintf( stderr,
             "-raw option invalid with more than one input file\n" );
    if( option->log )
    {
      writeLog( "-raw option invalid with more than one input file\n" );
    }
    option->errorFlag = 1;
  }

  if( option->edit == 0 && option->bottom )
  {
    fprintf( stderr,
             "-bottom option requires -view\n" );
    if( option->log )
    {
      writeLog( "-bottom option requires -view\n" );
    }
    option->errorFlag = 1;
  }
  if( option->edit == 0 && option->magnification )
  {
    fprintf( stderr,
             "-mag option requires -view\n" );
    if( option->log )
    {
      writeLog( "-mag option requires -view\n" );
    }
    option->errorFlag = 1;
  }
  if( option->edit == 0 && option->baseNumber )
  {
    fprintf( stderr,
             "-bn option requires -view\n" );
    if( option->log )
    {
      writeLog( "-bn option requires -view\n" );
    }
    option->errorFlag = 1;
  }
  if( option->edit == 0 && option->subSeq )
  {
    fprintf( stderr,
             "-loc option requires -view\n" );
    if( option->log )
    {
      writeLog( "-loc option requires -view\n" );
    }
    option->errorFlag = 1;
  }
  if( ( option->xmOption && option->call == 0 ) ||
      ( option->xmOption && option->edit == 0 ) )
  {
    fprintf( stderr, "-xm invalid without -view or with -nocall\n" );
    if( option->log )
    {
      writeLog( "-xm invalid without -view or with -nocall\n" );
    }
    option->errorFlag = 1;
  }
  if( !( option->scfVersion == 2 ||
         option->scfVersion == 3 ) )
  {
    fprintf( stderr, "-cv invalid argument value\n" );
    if( option->log )
    {
      writeLog( "-cv invalid argument value\n" );
    }
    option->errorFlag = 1;
  }
  if( !( option->scfPrecision == 1 ||
         option->scfPrecision == 2 ) )
  {
    fprintf( stderr, "-cp invalid argument value\n" );
    if( option->log )
    {
      writeLog( "-cp invalid argument value\n" );
    }
    option->errorFlag = 1;
  }
  if( option->trim == 1 &&
      option->call == 0 )
  {
    fprintf( stderr,
             "-trim option invalid with -nocall\n" );
    if( option->log )
    {
      writeLog( "-trim option invalid with -nocall\n" );
    }
    option->errorFlag = 1;
  }
  if( option->trim == 2 &&
      option->call == 0 )
  {
    fprintf( stderr,
             "-trim_alt option invalid with -nocall\n" );
    if( option->log )
    {
      writeLog( "-trim_alt option invalid with -nocall\n" );
    }
    option->errorFlag = 1;
  }
  if( option->trimFastaData == 1 && option->trim == 0 )
  {
    fprintf( stderr,
             "-trim_fasta option invalid without -trim or -trim_alt\n" );
    if( option->log )
    {
      writeLog( "-trim_fasta option invalid without -trim or -trim_alt\n" );
    }
    option->errorFlag = 1;
  }
  if( option->trimFastaData == 1 &&
      ( option->writeSeq == 0 && option->writeQual == 0 ) )
  {
    fprintf( stderr,
             "warning: -trim_fasta option set without -s, -sa, -sd, -q, -qa, or -qd\n" );
    if( option->log )
    {
      writeLog( "warning: -trim_fasta option set without -s, -sa, -sd, -q, -qa, or -qd\n" );
    }
  }
  if( option->trimFastaData == 1 && option->qualType != 0 )
  {
    fprintf( stderr,
             "-trim_fasta invalid with -qt xbap and -qt mix\n" );
    if( option->log )
    {
      writeLog( "-trim_fasta invalid with -qt xbap and -qt mix\n" );
    }
    option->errorFlag = 1;
  }
  if( option->trimSCFData == 1 && option->trim == 0 )
  {
    fprintf( stderr,
             "-trim_scf option invalid without -trim or -trim_alt\n" );
    if( option->log )
    {
      writeLog( "-trim_scf option invalid without -trim or -trim_alt\n" );
    }
    option->errorFlag = 1;
  }
  if( option->trimSCFData == 1 && option->writeScf == 0 )
  {
    fprintf( stderr,
             "warning: -trim_scf option set without -c or -cd\n" );
    if( option->log )
    {
      writeLog( "warning: -trim_scf set option without -c or -cd\n" );
    }
  }
  if( option->trimPHDData == 1 && option->trim == 0 )
  {
    fprintf( stderr,
             "-trim_phd option invalid without -trim or -trim_alt\n" );
    if( option->log )
    {
      writeLog( "-trim_phd option invalid without -trim or -trim_alt\n" );
    }
    option->errorFlag = 1;
  }
  if( option->trimPHDData == 1 && option->writePhd == 0 )
  {
    fprintf( stderr,
             "warning: -trim_phd option set without -p or -pd\n" );
    if( option->log )
    {
      writeLog( "warning: -trim_phd option set without -p or -pd\n" );
    }
  }
  if( option->dOption == 1 && option->call == 0 )
  {
    fprintf( stderr, "-d invalid with -nocall\n" );
    if( option->log )
    {
      writeLog( "-d invalid with -nocall\n" );
    }
    option->errorFlag = 1;
  }
  if( option->ddOption == 1 && option->call == 0 )
  {
    fprintf( stderr, "-dd invalid with -nocall\n" );
    if( option->log )
    {
      writeLog( "-dd invalid with -nocall\n" );
    }
    option->errorFlag = 1;
  }
  if( option->errorFlag )
  {
    return( ERROR );
  }
  else
  {
    return( OK );
  }
}
