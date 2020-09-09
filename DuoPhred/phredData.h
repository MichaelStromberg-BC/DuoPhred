/** phredData.h **/

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

#define FIL     	0
#define LCL     	1

typedef struct
{
  int    alloc;
  int    numBase;
  FLOAT  A_scale_factor;
  FLOAT  C_scale_factor;
  FLOAT  G_scale_factor;
  FLOAT  T_scale_factor;
  int    *called_nuc;
  int    *called_loc;
  FLOAT  *called_area;
  FLOAT  *called_relative_area;
  int    *uncalled_nuc;
  int    *uncalled_loc;
  FLOAT  *uncalled_area;
  FLOAT  *uncalled_relative_area;
  FLOAT  *A_tr_val;
  FLOAT  *C_tr_val;
  FLOAT  *G_tr_val;
  FLOAT  *T_tr_val;
} PolyData;

typedef struct
{
  char   fileName[PHRED_PATH_MAX];
  int    fileType;

  int    chemType;
  int    dyeType;
  int    machineType;

  int    compressSplitFlag;
  int    primerQVFlag;

  int    numPoint;
  FLOAT *trace[4];

  int    numBase[2];
  char  *base[2];
  char  *baseUncalled[2];
  int   *baseLoc[2];
  int   *baseLocUncalled[2];
  int   *baseQual[2];
  int   *baseQualUncalled[2];
  int    numDsc;
  int   *baseDsc;
  int    begAlign;
  int    endAlign;
  FLOAT *qualIndex;

  int    dataSet[2];

  FLOAT  maxTraceValue;
  FLOAT  minTraceValue;

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
 
  int    leftTrimPoint;
  int    rghtTrimPoint;
  int    lftPhredTrim;
  int    rhtPhredTrim;
  FLOAT  trimSetValue;
  int    qualType;

  PolyData polyData;

} PhredData;


