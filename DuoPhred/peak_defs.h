/** peak_defs.h **/


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
#include <math.h>
#include <ctype.h>

#define MAX_NUM_SPLITS 4 /* maximum number of peaks into which a single
			    observed peak can be split */
#define MAX_RIGHT_SHIFT 2 /* maximum right shift allowed (as number of
			     predicted peaks */
#define MAX_LEFT_SHIFT 4
#define SPLIT_PENALTY .8 /* penalty for splitting peak. 0 corresponds to no
			  penalty; 1 to complete penalty (0 score for
			  additional peaks) */
#define LEFT_SHIFT_PENALTY .1
#define RIGHT_SHIFT_PENALTY .3
#define MAX_SHIFT_CHANGE .7
#define LEFT_SPLIT_SHIFT_PENALTY .1
#define RIGHT_SPLIT_SHIFT_PENALTY .3
#define MAX_LEFT_DISPLACEMENT 2.1
#define MAX_RIGHT_DISPLACEMENT .5
#define MIN_PEAK_RATIO .10
#define RIGHT_DISPLACEMENT_WINDOW .5 
#define RIGHT_EXPANDED_WINDOW .5 
#define LEFT_DISPLACEMENT_WINDOW .5
#define LEFT_EXPANDED_WINDOW .6
#define MAX_COMPRESSION_LENGTH 5
#define MIN_RELATIVE_AREA .2
#define MIN_SPLITTABLE_AREA 1.6 /* minimum relative area necessary to permit
				   positive score splitting of peak */
				   
typedef struct peak {
  FLOAT pred_location, pred_period, proportion_fitted;
  FLOAT total_signal;
  FLOAT last10; /* average strength of preceding 10 peaks */
  char nuc; /* the called nucleotide */
  struct observed_peak *obs_peak, *best_obs_peak, *best_uncalled_peak;
  int next_no_observed; 
  char fixed; /* indicates that peak was assigned in the loop which assigns
		evenly spaced consecutive observed peaks to peaks, in
		fit_peaks */
  struct peak *next, *prev; /* doubly linked list */
} Peak;

typedef struct observed_peak {
  FLOAT area, relative_area, relative_nuc_area;
  char nuc; /* 0,1,2,3 */
  int split[MAX_NUM_SPLITS + 1][MAX_NUM_SPLITS + 1];
  int first_peak_index, last_peak_index;
  char dp[MAX_RIGHT_SHIFT + MAX_LEFT_SHIFT + 1];  /* dynamic programming vector: */
  FLOAT scores[MAX_RIGHT_SHIFT + MAX_LEFT_SHIFT + 1];
  FLOAT shift;
  int i_left;
  int i_rite;
  int i_maxx;
  int type;
  Peak *peak[MAX_NUM_SPLITS + 1];
  struct observed_peak *next, *prev;
} Observed_peak;

typedef struct cos_sin_array {
  int length;
  double period;
  double *cos, *sin;
  struct cos_sin_array *next;
} Cos_sin_array;

typedef struct
{
  int   loc;
  int   type;
  int   nuc;
  int   basMaxSpc;
  int   locUncalled;
  int   nucUncalled;
  int   flgUncalled;
  FLOAT maxCalled;
  FLOAT maxUncalled;
  FLOAT omaxCalled;
  FLOAT omaxUncalled;
  FLOAT maxRatio3;
  FLOAT maxRatio7;
  FLOAT maxRatio3Uncalled;
  FLOAT maxRatio7Uncalled;
  FLOAT relativeArea;
  FLOAT space;
  FLOAT spcRatio;
  FLOAT maxDownF1;
  FLOAT maxDown;
  FLOAT basAreaVar;
  FLOAT maxAreaVar;
  FLOAT compPar;
  FLOAT compMaxArea;
  FLOAT resPar7;
  Peak  *peak;
} BaseQual;

typedef struct
{
  int base;
  int point;
  int quality;
  int type;
  int width;
  int shoulder;
  int truncated;
  FLOAT stdev;
  FLOAT space;
  FLOAT spcRatio;
} LocPeak;

#ifdef ANSI_C

int fit_main( int tr_length, FLOAT **tr_vals,
              PhredData *phredData, int *status );
Peak *fit_sine( FLOAT **tr_vals, FLOAT *tot_vals, int tr_length,
                FLOAT *scale_traces, Observed_peak **pfirst_obs_peak,
                int compressSplitFlag, int chemType, int *status );

#else

int fit_main();
Peak *fit_sine();

#endif
