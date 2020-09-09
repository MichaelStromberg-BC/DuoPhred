/** phred.c **/

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

#ifdef X86_GCC_LINUX
#include <fpu_control.h>
#endif

#ifdef ANSI_C
int main( int argc, char **argv )
#else
int main( argc, argv )
int argc;
char **argv;
#endif
{
  Option option;

#ifdef X86_GCC_LINUX
  /*
  ** Set the x86 FPU control word to force double
  ** precision rounding rather than `extended'
  ** precision rounding. This causes phred base
  ** calls and quality values on x86 GCC-Linux
  ** (tested on RedHat Linux) machines to be
  ** identical to those on IEEE conforming UNIX
  ** machines.
  */
  fpu_control_t cw;

  cw = ( _FPU_DEFAULT & ~_FPU_EXTENDED ) | _FPU_DOUBLE;

  _FPU_SETCW( cw );
#endif

  /*
  ** provide access to option structure
  */
  setOption( &option );

  /*
  ** initialize parameters
  */
  initParam( &option, argc );

  /*
  ** read .phredrc
  */
  readRC( &option );

  /*
  ** get environment variables
  */
  readEnvVar( &option );

  /*
  ** get command line parameters
  */
  readParam( argc, argv, &option );

  /*
  ** open log file if requested
  */
  if( option.log )
  {
    if( openLog( ) == ERROR )
    {
      return( ERROR );
    }
    logVersion();
    logParam( argc, argv );
  }

  /*
  ** check parameter consistency
  */
  checkParam( &option );

  /*
  ** quit if an error occurred
  */
  if( option.errorFlag )
  {
    if( option.log )
    {
      writeLog( "\n*** Bad parameters...quitting ***\n" );
      closeLog();
    }
    freeParam( &option );
    return( ERROR );
  }

  /*
  ** Read phred parameter file, if necessary.
  */
  if( option.parFileName[0] != '\0' )
  {
    if( readParamFile( &option ) == ERROR )
    {
      fprintf( stderr, "  warning: processing without phred parameters\n" );
    }
  }

  /*
  ** Open files for concatenated output.
  */
  if( option.saOption )
  {
    if( ( option.seqFP = fopen( option.writeSeqName, "w+" ) ) == NULL )
    {
      fprintf( stderr, "unable to open file %s\n", option.writeSeqName );
      if( option.log )
      {
        writeLog( "phred: error: unable to open sequence FASTA file\n" );
      }
    }
  }
  if( option.qaOption )
  {
    if( ( option.qualFP = fopen( option.writeQualName, "w+" ) ) == NULL )
    {
      fprintf( stderr, "unable to open file %s\n", option.writeQualName );
      if( option.log )
      {
        writeLog( "phred: error: unable to open quality FASTA file\n" );
      }
    }
  }

  /*
  ** choose a PHRED
  */
  if( option.edit )
  {
    if( viewPhred( argc, argv, &option )  == ERROR )
    {
      if( option.log )
      {
        writeLog( "phred: error: bad status: viewPhred\n" );
        closeLog();
      }
      freeParam( &option );
      return( ERROR );
    }
  }
  else
  {
    if( autoPhred( &option ) == ERROR )
    {
      if( option.log )
      {
        writeLog( "phred: error: bad status: autoPhred\n" );
        closeLog();
      }
      freeParam( &option );
      return( ERROR );
    }
  }

  /*
  ** close log file if opened
  */
  if( option.log )
  {
    closeLog();
  }

  /*
  ** Close appended output files.
  */
  if( option.saOption )
  {
    fclose( option.seqFP );
  }
  if( option.qaOption )
  {
    fclose( option.qualFP );
  }

  /*
  ** free memory for parameters
  */
  freeParam( &option );

  return( OK );
}


