/** readESD.h **/

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

#include "chromatData.h"

typedef unsigned int4           DWORD;
typedef unsigned char           TCHAR;
typedef char                    BYTE;
typedef int4                    BOOL;
typedef unsigned int            UINT;
typedef int                     ScanRate;
typedef int4                    l_time_t;

#define INSTRUMENT_NAME_LEN     ((10*2+2)/sizeof(TCHAR))
#define SAMPLE_NAME_LEN         (24 / sizeof(TCHAR))
#define SPLITTER_NAME_LEN       (30 / sizeof(TCHAR))
#define FILTER_NAME_LEN         (30 / sizeof(TCHAR))
#define OPERATOR_NAME_LEN       (30 / sizeof(TCHAR))
#define DYE_NAME_LEN            (10 / sizeof(TCHAR))
#define PAR_SET_NAME_LEN        (30 / sizeof(TCHAR))
#define PLATE_NAME_LEN          (18 / sizeof(TCHAR))
#define RUN_ID_LEN              (6  / sizeof(TCHAR))
#define CHEM_NAME_LEN           (32 / sizeof(TCHAR))

#define ORDER_A                 0
#define ORDER_C                 1
#define ORDER_G                 2
#define ORDER_T                 3

#define ES_FOUR_CHAN            4

#define ESD_STR_MAX		256

typedef struct
{
  int4 capillaryCurrent;
  float channel1Intensity;
  float channel2Intensity;
  float channel3Intensity;
  float channel4Intensity;
} ProcessedScan;


typedef struct
{
  TCHAR         m_tcSampleName[SAMPLE_NAME_LEN];
  TCHAR         m_tcOperator[OPERATOR_NAME_LEN];

  float         m_fRunVoltage;
  float         m_fRunTime;
  float         m_fInjectionTime;
  float         m_fInjectionVoltage;
  float         m_fTemperature;

  int           m_nLazerIndex;
  float         m_fPMTVoltage1;
  float         m_fPMTVoltage2;
  TCHAR         m_csBeamSplitterA[SPLITTER_NAME_LEN];
  TCHAR         m_csBeamSplitterB[SPLITTER_NAME_LEN];
  TCHAR         m_csFilter1[FILTER_NAME_LEN];
  TCHAR         m_csFilter2[FILTER_NAME_LEN];
  TCHAR         m_csFilter3[FILTER_NAME_LEN];
  TCHAR         m_csFilter4[FILTER_NAME_LEN];

  float         m_fFillTime;
  float         m_fFlushTime1;
  float         m_fFlushTime2;
  float         m_fFlushTime3;
  float         m_fPrerunVoltage;
  float         m_fRelaxationTime;
  float         m_fPrerunTime;
  float         m_fSpecSepMx[ES_FOUR_CHAN][ES_FOUR_CHAN];

  TCHAR         m_csBaseOrder[ES_FOUR_CHAN][DYE_NAME_LEN];
  TCHAR         m_InstID[INSTRUMENT_NAME_LEN];

  TCHAR         m_csParamSetName[PAR_SET_NAME_LEN];

  DWORD         m_dMagicNumber;                 /* 826561349 */
  DWORD         m_dExtHeaderOffset;
  DWORD         m_dExtHeaderLen;

  float         m_fPreinjectionTime;
  float         m_fPreinjectionVoltage;
  float         m_fLowPressureTime;
  float         m_fUserInputTime;
  float         m_fSleepTime;
  float         m_fSleepTemperature;

  BYTE          m_byReserved[82];

  TCHAR         m_csRunID[RUN_ID_LEN];
  TCHAR         m_csPlateID[PLATE_NAME_LEN];
  DWORD         m_dwPrevRunID;

  BOOL          m_bDWORDFile;

  BOOL          m_bSaturated;
  float         m_fVersion;
  l_time_t      m_tRunStartTime;
  l_time_t      m_tRunStopTime;

  DWORD         m_dwNumberOfLines;

  ScanRate      m_nScanRate;

  BYTE          m_byFileType;

  BYTE          m_byNumberOfCapillaries;

  BYTE          m_byCapUsage;
  BYTE          m_byDummy;

  UINT          m_uiSecurityCode;

} ESDTrailer;

typedef struct
{
  char primerName[ESD_STR_MAX];
  char sampleName[ESD_STR_MAX];
  char plateID[ESD_STR_MAX];
  char comment[ESD_STR_MAX];
  char machineName[ESD_STR_MAX];
  char wellID[ESD_STR_MAX];
  char baseMap[4];
  int  numBase;
  char *base;
  int  *baseLoc;
  int  baseLocOffset;
} ESDExtHeader;


#ifdef ANSI_C
ChromatData *readESD( char *filename, int version, int *status );
#else
ChromatData *readESD();
#endif
