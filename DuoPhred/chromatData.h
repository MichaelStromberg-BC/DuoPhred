/** chromatData.h **/

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
** Data from chromat file.
*/

#ifndef CHROMATDATA_DEFINED
#define CHROMATDATA_DEFINED


#ifndef FLOAT
#define FLOAT		double
#endif

/*
** This should be only in typeDef.h but is
** repeated here for David's convenience.
** The same applies to DOUBLE above.
*/
#ifdef PHRED_PATH_MAX
#undef PHRED_PATH_MAX
#endif

#define PHRED_PATH_MAX	4096

#ifndef SCFFormat
#define SCFFormat	1
#endif

#ifndef ABIFormat
#define ABIFormat	2
#endif

#ifndef MD1Format
#define MD1Format       3
#endif

#ifndef MD2Format
#define MD2Format       4
#endif

typedef struct
{
  char   fileName[PHRED_PATH_MAX];
  int    fileType;

  int    numPoint;
  FLOAT *trace[4];
  FLOAT  maxTraceValue;
  FLOAT  minTraceValue;

  int    numBase;
  char  *base;
  char	*baseUncalled;
  int   *baseLoc;
  int	*baseLocUncalled;
  int   *baseQual;
  int	*baseQualUncalled;

  int    laneNumber;
  int    primerLoc;
  int    signalStrength[4];
  FLOAT  avgSpacing;
  char   comment[PHRED_MAX_STRING_LEN];
  char   convProg[PHRED_MAX_STRING_LEN];
  char   gelName[PHRED_MAX_STRING_LEN];
  char   machineName[PHRED_MAX_STRING_LEN];
  char   primerID[PHRED_MAX_STRING_LEN];
  char   processing[PHRED_MAX_STRING_LEN];
  char   reTracker[PHRED_MAX_STRING_LEN];
  char   sampleName[PHRED_MAX_STRING_LEN];
  char   source[PHRED_MAX_STRING_LEN];

  char   thumbPrint[10];

} ChromatData;

#ifdef ANSI_C
ChromatData *allocChromatData( int numPoint, int numBase );
int freeChromatData( ChromatData *chromatData );
#else
ChromatData *allocChromatData();
int freeChromatData();
#endif

#endif
