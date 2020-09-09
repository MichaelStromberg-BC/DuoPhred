/** setQual.c **/

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

typedef struct
{
  int quality;
  FLOAT maxDown;
  FLOAT ampRatio3;
  FLOAT ampRatio7;
  FLOAT spcRatio;
} QualValueFourPar;

typedef struct
{
  int quality;
  FLOAT maxDown;
  FLOAT ampRatio3;
  FLOAT ampRatio7;
  FLOAT spcRatio;
  FLOAT compPar;
} QualValueFivePar;


/*
** Include quality value lookup tables.
*/
#include "qualTable1.h"
#include "qualTable2.h"
#include "qualTable3.h"
#include "qualTable4.h"
#include "qualTable5c.h"
#include "qualTable6c.h"

#ifdef ANSI_C
int setQual( int tr_length, FLOAT **tr_vals, FLOAT *tot_vals,
             int num_base, Peak *first_peak, int *quality, int *qualityUncalled, PhredData *phredData,
             int qualityValueCeilingOption, int qualityValueCeiling )
#else
int setQual( tr_length, tr_vals, tot_vals, num_base, first_peak, quality,
             phredData, qualityValueCeilingOption, qualityValueCeiling )
int tr_length;
FLOAT **tr_vals;
FLOAT *tot_vals;
int num_base;
Peak *first_peak;
int *quality;
int qualityValueCeilingOption;
int qualityValueCeiling;
PhredData *phredData;
#endif
{
  int i, j, z;
  int numBase;
  BaseQual *baseQual;
  Option *option;
	float mpsMaxRatio;
	//FLOAT mpsNextSpace;
	FILE* OUT;

	int mpsUncalledNuc;
	int mpsUncalledLoc;
	int mpsCalledNuc;
	int mpsCalledLoc;
	int uncalledQuality;
	
  option = getOption();

  if( num_base < 10 )
  {
    for( i = 0; i < num_base; ++i )
    {
      quality[i] = 0;
    }
    return( OK );
  }

  /*
  ** Initialize the base quality structure array.
  */
  baseQual = initBaseQual( first_peak, &numBase );

  if( num_base != numBase )
  {
    printf( "setQual: inconsistent number of bases\n" );
  }

  /*
  ** Set the maxDown elements of baseQual.
  */
  maxDownBaseQual( tr_length, tr_vals, tot_vals, numBase, baseQual );

  /*
  ** Set the space ratio elements of baseQual.
  */
  spcRatioBaseQual( numBase, baseQual );

  /*
  ** Set the maxRatio elements of baseQual.
  */
  if( phredData->chemType == TERMINATOR_CHEM &&
      phredData->dyeType  == RHODAMINE_DYE )
  {
    maxRatioBaseQual( tr_length, tr_vals, phredData->trace, 1, numBase, baseQual );
  }
  else
  {
  	//printf("MPS: Calling maxRatioBaseQual. tr_length: %d, numBase: %d\n",tr_length,numBase);
    maxRatioBaseQual( tr_length, tr_vals, phredData->trace, 0, numBase, baseQual );
  }
  
  if( phredData->primerQVFlag == 1 && phredData->chemType == PRIMER_CHEM )
  {
    /*
    ** Set the compPar elements of baseQual.
    */
    compressParameter( numBase, baseQual );
  }

  /*
  ** Calculate the base qualities.
  */
  if( phredData->machineType == ABI_373_377 ||
      phredData->machineType == ABI_3700    ||
      phredData->machineType == UNKNOWN_MACHINE )
  {

    if( option->verboseOption && option->verboseLevel >= 1 )
    {
      fprintf( stderr,
               "setQual: ABI/UNKNOWN quality values: " );
    }

    if( phredData->primerQVFlag == 1 && phredData->chemType == PRIMER_CHEM )
    {

       if( option->verboseOption && option->verboseLevel >= 1 )
       {
         fprintf( stderr,
                  "five parameter for dye primer\n" );
       }

      /*
      ** Primer chemistry specified.
      */
      for( i = 0; i < numBase; ++i )
      {
        quality[i] = 0;
        if( baseQual[i].nuc == 4 )
        {
          continue;
        }
        for( j = 0; j < NMQL2; ++j )
        {
          if( baseQual[i].maxDown    <= qualLim2[j].maxDown   &&
              baseQual[i].maxRatio3  <= qualLim2[j].ampRatio3 &&
              baseQual[i].maxRatio7  <= qualLim2[j].ampRatio7 &&
              baseQual[i].spcRatio   <= qualLim2[j].spcRatio  &&
              baseQual[i].compPar    <= qualLim2[j].compPar )
          {
            quality[i] = qualLim2[j].quality;
          }
        }
        if( qualityValueCeilingOption &&
            quality[i] > qualityValueCeiling )
        {
          quality[i] = qualityValueCeiling;
        }
      }
    }
    else
    {

			//
			// MPS: This is what we're using
			//

      if( option->verboseOption && option->verboseLevel >= 1 )
      {
        fprintf( stderr,
                 "four parameter for dye terminator and unknown\n" );
      }

      /*
      ** Four parameter quality values.
      */
      for( i = 0; i < numBase; ++i )
      {
        quality[i]         = 0;
		qualityUncalled[i] = 0;
        
        //
        // handle the called base
        //
        if(baseQual[i].nuc == 4) continue;
       
        for( j = 0; j < NMQL1; ++j ) {
          if( baseQual[i].maxDown    <= qualLim1[j].maxDown   &&
              baseQual[i].maxRatio3  <= qualLim1[j].ampRatio3 &&
              baseQual[i].maxRatio7  <= qualLim1[j].ampRatio7 &&
              baseQual[i].spcRatio   <= qualLim1[j].spcRatio  ) {
            quality[i] = qualLim1[j].quality;
          }
        }
        
        //
        // handle the uncalled base
        //
  
        if(baseQual[i].nucUncalled < 4) {
  				        	
        	for( j = 0; j < NMQL1; ++j ) {
	          if( baseQual[i].maxDown            <= qualLim1[j].maxDown   &&
	              baseQual[i].maxRatio3Uncalled  <= qualLim1[j].ampRatio3 &&
	              baseQual[i].maxRatio7Uncalled  <= qualLim1[j].ampRatio7 &&
	              baseQual[i].spcRatio           <= qualLim1[j].spcRatio  ) {
	            qualityUncalled[i] = qualLim1[j].quality;
	          }
       		}
     
     			//printf("           %d, uncalled [%c] bq: %d\n",i,"ACGTN"[baseQual[i].nucUncalled],uncalledQuality);
        }
        
        
        if( qualityValueCeilingOption &&
            quality[i] > qualityValueCeiling )
        {
          quality[i] = qualityValueCeiling;
        }
	  }

		//
		// MPS: Save our XML file now
		//
      if( option->xdOption == 1) writeXmlData(num_base,quality,qualityUncalled,baseQual,phredData);
    }
  }
  else
  if( phredData->machineType == MOLDYN_MEGABACE )
  {

    if( option->verboseOption && option->verboseLevel >= 1 )
    {
      fprintf( stderr,
               "setQual: MegaBACE quality values: " );
    }

    if( phredData->primerQVFlag == 1 && phredData->chemType == PRIMER_CHEM )
    {

      if( option->verboseOption && option->verboseLevel >= 1 )
      {
        fprintf( stderr,
                 "five parameter for dye primer\n" );
      }

      /*
      ** Primer chemistry specified.
      */
      for( i = 0; i < numBase; ++i )
      {
        quality[i] = 0;
        if( baseQual[i].nuc == 4 )
        {
          continue;
        }
        for( j = 0; j < NMQL4; ++j )
        {
          if( baseQual[i].maxDown    <= qualLim4[j].maxDown   &&
              baseQual[i].maxRatio3  <= qualLim4[j].ampRatio3 &&
              baseQual[i].maxRatio7  <= qualLim4[j].ampRatio7 &&
              baseQual[i].spcRatio   <= qualLim4[j].spcRatio  &&
              baseQual[i].compPar    <= qualLim4[j].compPar )
          {
            quality[i] = qualLim4[j].quality;
          }
        }
        if( qualityValueCeilingOption &&
            quality[i] > qualityValueCeiling )
        {
          quality[i] = qualityValueCeiling;
        }
      }
    }
    else
    {

      if( option->verboseOption && option->verboseLevel >= 1 )
      {
        fprintf( stderr,
                 "four parameter for dye terminator and unknown\n" );
      }

      /*
      ** Four parameter quality values.
      */
      for( i = 0; i < numBase; ++i )
      {
        quality[i] = 0;
        if( baseQual[i].nuc == 4 )
        {
          continue;
        }
        for( j = 0; j < NMQL3; ++j )
        {
          if( baseQual[i].maxDown    <= qualLim3[j].maxDown   &&
              baseQual[i].maxRatio3  <= qualLim3[j].ampRatio3 &&
              baseQual[i].maxRatio7  <= qualLim3[j].ampRatio7 &&
              baseQual[i].spcRatio   <= qualLim3[j].spcRatio  )
          {
            quality[i] = qualLim3[j].quality;
          }
        }
        if( qualityValueCeilingOption &&
            quality[i] > qualityValueCeiling )
        {
          quality[i] = qualityValueCeiling;
        }
      }
    }
  }
  else
  if( phredData->machineType == LI_COR_4000 )
  {

    if( option->verboseOption && option->verboseLevel >= 1 )
    {
      fprintf( stderr,
               "setQual: LI-COR quality values: " );
    }

    if( phredData->primerQVFlag == 1 && phredData->chemType == PRIMER_CHEM )
    {
      if( option->verboseOption && option->verboseLevel >= 1 )
      {
        fprintf( stderr,
                 "five parameter for dye primer\n" );
      }

      /*
      ** Primer chemistry specified.
      */
      for( i = 0; i < numBase; ++i )
      {
        quality[i] = 0;
        if( baseQual[i].nuc == 4 )
        {
          continue;
        }
        for( j = 0; j < NMQL6; ++j )
        {
          if( baseQual[i].maxDown    <= qualLim6[j].maxDown   &&
              baseQual[i].maxRatio3  <= qualLim6[j].ampRatio3 &&
              baseQual[i].maxRatio7  <= qualLim6[j].ampRatio7 &&
              baseQual[i].spcRatio   <= qualLim6[j].spcRatio  &&
              baseQual[i].compPar    <= qualLim6[j].compPar )
          {
            quality[i] = qualLim6[j].quality;
          }
        }
        if( qualityValueCeilingOption &&
            quality[i] > qualityValueCeiling )
        {
          quality[i] = qualityValueCeiling;
        }
      }
    }
    else
    {

      if( option->verboseOption && option->verboseLevel >= 1 )
      {
        fprintf( stderr,
                 "four parameter for dye terminator and unknown\n" );
      }

      /*
      ** Four parameter quality values.
      */
      for( i = 0; i < numBase; ++i )
      {
        quality[i] = 0;
        if( baseQual[i].nuc == 4 )
        {
          continue;
        }
        for( j = 0; j < NMQL5; ++j )
        {
          if( baseQual[i].maxDown    <= qualLim5[j].maxDown   &&
              baseQual[i].maxRatio3  <= qualLim5[j].ampRatio3 &&
              baseQual[i].maxRatio7  <= qualLim5[j].ampRatio7 &&
              baseQual[i].spcRatio   <= qualLim5[j].spcRatio  )
          {
            quality[i] = qualLim5[j].quality;
          }
        }
        if( qualityValueCeilingOption &&
            quality[i] > qualityValueCeiling )
        {
          quality[i] = qualityValueCeiling;
        }
      }
    }
  }
  else
  {
    fprintf( stderr,
             "setQual: unrecognized machine type specification\n" );
    ourFree( (char *)baseQual );
    return( ERROR );
  }




  /*
  ** Free memory.
  */
  ourFree( (char *)baseQual );

  return( OK );
}


