/** readParamFile.c **/

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
int readParamFile( Option *option )
#else
int readParamFile( option )
Option *option;
#endif
{
  int i;
  int iset;
  int istate;
  int iline;
  int nline;
  int numByte;
  int numChem;
  char *cptr;
  char line[1024];
  char set_id[1024];
  ChemistryList *chemList;
  FILE *fp;

  /*
  ** Open parameter file.
  */
  fp = fopen( option->parFileName, "r" );
  if( fp == NULL )
  {
    fprintf( stderr,
             "readParamFile: warning: unable to open file %s\n",
             option->parFileName );
    return( ERROR );
  }

  /*
  ** Count lines.
  */
  nline = 0;
  while( fgets( line, 1024, fp ) != NULL )
  {
    /*
    ** Skip leading whitespace.
    */
    i = 0; 
    while( line[i] == ' ' ||
           line[i] == '\t' )
    {
      ++i;
    }

    cptr = &(line[i]);

    /*
    ** Skip blank lines.
    */
    if( *cptr == '\n' )
    {
      continue;
    }

    /*
    ** Skip comment lines.
    */
    if( *cptr == '#' )
    {
      continue;
    }

    ++nline;
  }

  /*
  ** Read parameter entries.
  */
  istate = 0; 
  iset   = 0;
  iline  = 0;
  rewind( fp );
  while( fgets( line, 1024, fp ) != NULL )
  {
    ++iline;

    /*
    ** Skip leading whitespace.
    */
    i = 0;
    while( line[i] == ' ' ||
           line[i] == '\t' )
    {
      ++i;
    }

    cptr = &(line[i]);

    /*
    ** Skip blank lines.
    */
    if( *cptr == '\n' )
    {
      continue;
    }

    /*
    ** Skip comment lines.
    */
    if( *cptr == '#' )
    {
      continue;
    }

    if( istate == 1 )
    {
      if( strncmp( cptr, "end", 3 ) == 0 )
      {
        if( iset == 1 )
        {
          numByte = numChem * sizeof( ChemistryList );
          option->parFileData.chemList = (ChemistryList *)malloc( numByte );
          if( option->parFileData.chemList == NULL )
          {
            fprintf( stderr, "readParamFile: error: unable to allocate memory\n" );
            fclose( fp );
            return( ERROR );
          }
          memcpy( option->parFileData.chemList, chemList, numByte );

          /*
          ** Free chemList memory but not the `name' strings.
          */
          free( chemList );
        }
        option->parFileData.numChem = numChem;
        istate = 0;
        iset   = 0;
      }
      else
      if( iset == 1 )
      {
        /*
        ** Parse out chemistry list lines.
        ** The primer name must be enclosed in double quotes.
        */
        if( *cptr != '"' )
        {
          fprintf( stderr, "readParamFile: error: line %d missing `\"'\n", iline );
          freeChemList( nline, chemList );
          fclose( fp );
          return( ERROR );
        }

        ++cptr;

        i = 0;
        while( cptr[i] != '"' )
        {
          if( cptr[i] == '\0' )
          {
            fprintf( stderr, "readParamFile: error: line %d missing `\"'\n", iline );
            freeChemList( nline, chemList );
            fclose( fp );
            return( ERROR );
          }
          ++i;
        }
        cptr[i] = '\0';

        numByte = ( strlen( cptr ) + 1 ) * sizeof( char );
        chemList[numChem].name = (char *)malloc( numByte );
        if( chemList[numChem].name == NULL )
        {
          fprintf( stderr, "readParamFile: error: unable to allocate memory\n" );
          freeChemList( nline, chemList );
          fclose( fp );
          return( ERROR );
        }
        strcpy( chemList[numChem].name, cptr );

        /*
        ** Parse out chemistry type.
        */
        cptr = cptr + i + 1;

        i = 0;
        while( cptr[i] == ' ' ||
               cptr[i] == '\t' )
        {
          ++i;
        }

        cptr = cptr + i;
        i = 0;
        while( cptr[i] != ' ' &&
               cptr[i] != '\t' &&
               cptr[i] != '\n' &&
               cptr[i] != '\0' )
        {
          ++i;
        }
        if( i == 0 ||
            cptr[i] == '\n' ||
            cptr[i] == '\0' )
        {
          fprintf( stderr,
                   "readParamFile: error: line %d missing field\n",
                   iline );
          freeChemList( nline, chemList );
          fclose( fp );
          return( ERROR );
                   
        }
        cptr[i] = '\0';

        if( strncmp( cptr, "primer", 6 ) == 0 )
        {
          chemList[numChem].chemType = PRIMER_CHEM;
        }
        else
        if( strncmp( cptr, "terminator", 9 ) == 0 )
        {
          chemList[numChem].chemType = TERMINATOR_CHEM;

        }
        else
        if( strncmp( cptr, "unknown", 7 ) == 0 )
        {
          chemList[numChem].chemType = UNKNOWN_CHEM;
        }
        else
        {
          fprintf( stderr,
                   "readParamFile: warning: line %d unknown chemistry: %s\n",
                   iline, cptr );
          chemList[numChem].chemType = UNKNOWN_CHEM;
        }

        /*
        ** Parse out dye type.
        */
        cptr = cptr + i + 1;

        i = 0;
        while( cptr[i] == ' ' ||
               cptr[i] == '\t' )
        {
          ++i;
        }

        cptr = cptr + i;
        i = 0;
        while( cptr[i] != ' ' &&
               cptr[i] != '\t' &&
               cptr[i] != '\n' &&
               cptr[i] != '\0' )
        {
          ++i;
        }
        if( i == 0 ||
            cptr[i] == '\n' ||
            cptr[i] == '\0' )
        {
          fprintf( stderr,
                   "readParamFile: error: line %d missing field\n",
                   iline );
          freeChemList( nline, chemList );
          fclose( fp );
          return( ERROR );
                   
        }
        cptr[i] = '\0';

        if( strncmp( cptr, "orig-dye",  8 ) == 0 ||
            strncmp( cptr, "rhodamine", 9 ) == 0 )
        {
          chemList[numChem].dyeType = RHODAMINE_DYE;
        }
        else
        if( strncmp( cptr, "big-dye", 7 ) == 0 )
        {
          chemList[numChem].dyeType = BIG_DYE_DYE;
        }
        else
        if( strncmp( cptr, "d-rhodamine", 11 ) == 0 )
        {
          chemList[numChem].dyeType = D_RHODAMINE_DYE;
        }
        else
        if( strncmp( cptr, "energy-transfer", 15 ) == 0 )
        {
          chemList[numChem].dyeType = ENERGY_TRANSFER_DYE;
        }
        else
        if( strncmp( cptr, "bodipy", 6 ) == 0 )
        {
          chemList[numChem].dyeType = BODIPY_DYE;
        }
        else
        if( strncmp( cptr, "unknown", 7 ) == 0 )
        {
          chemList[numChem].dyeType = UNKNOWN_DYE;
        }
        else
        {
          fprintf( stderr,
                   "readParamFile: error: %s: unknown dye type: %s\n",
                   option->parFileName,
                   cptr );
          chemList[numChem].dyeType = UNKNOWN_DYE;
        }

        /*
        ** Parse out machine type.
        */
        cptr = cptr + i + 1;

        i = 0;
        while( cptr[i] == ' ' ||
               cptr[i] == '\t' )
        {
          ++i;
        }

        cptr = cptr + i;
        i = 0;
        while( cptr[i] != ' ' &&
               cptr[i] != '\t' &&
               cptr[i] != '\n' &&
               cptr[i] != '\0' )
        {
          ++i;
        }
        if( i == 0 )
        {
          fprintf( stderr,
                   "readParamFile: error: line %d missing field\n",
                   iline );
          freeChemList( nline, chemList );
          fclose( fp );
          return( ERROR );
        }
        cptr[i] = '\0';


        if( strncmp( cptr, "ABI_373_377", 11 ) == 0 )
        {
          chemList[numChem].machineType = ABI_373_377;
        }
        else
        if( strncmp( cptr, "MolDyn_MegaBACE", 15 ) == 0 )
        {
          chemList[numChem].machineType = MOLDYN_MEGABACE;
        }
        else
        if( strncmp( cptr, "ABI_3700", 8 ) == 0 )
        {
          chemList[numChem].machineType = ABI_3700;
        }
        else
        if( strncmp( cptr, "LI-COR_4000", 11 ) == 0 )
        {
          chemList[numChem].machineType = LI_COR_4000;
        }
        else
        {
          fprintf( stderr,
                   "readParamFile: error: %s: unknown machine type: %s\n",
                   option->parFileName,
                   cptr );
          chemList[numChem].machineType = UNKNOWN_MACHINE;
        }

        ++numChem;
      }
    }
    else
    if( istate == 0 )
    {
      if( strncmp( cptr, "begin", 5 ) == 0 )
      {
        if( sscanf( &(cptr[6]), "%s", set_id ) != 1 )
        {
          fprintf( stderr,
                   "readParamFile: warning: %d: `begin' missing set type\n",
                   iline );
          istate = 1;
          iset   = 0;
          continue;
        }

        if( strncmp( set_id, "chem_list", 9 ) == 0 )
        {
          numByte = nline * sizeof( ChemistryList );
          chemList = (ChemistryList *)malloc( numByte );
          if( chemList == NULL )
          {
            fprintf( stderr, "readParamFile: error: unable to allocate memory\n" );
            fclose( fp );
            return( ERROR );
          }
          for( i = 0; i < nline; ++i )
          {
            chemList[i].name = NULL;
          }
          iset = 1;
          numChem = 0;
        }
        else
        {
          fprintf( stderr,
                   "readParamFile: warning: line %d: `begin' unknown set type %s\n",
                   iline,
                   set_id );
          iset = 0;
        }
        istate = 1;
      }
    }
  }

  if( istate == 1 )
  {
    if( iset == 1 )
    {
      freeChemList( nline, chemList );
    }
    fprintf( stderr, "readParamFile: error: missing end for %s\n", set_id );
    fclose( fp );
    return( ERROR );
  }

  /*
  ** Note a successful file read.
  */
  option->parFileData.successfulRead = 1;

  /*
  ** Close file.
  */
  fclose( fp );

  return( OK );
}

