/** typeDef.h **/

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
** Take care of definitions that are repeated in chromatData.h for
** David's convenience.
*/

#define OK	0
#define	ERROR	-1

/*
** writeSeq.c uses these format codes to indicate
** the input file format in the header of the
** output sequence files.
*/
#define EEKFormat 	-1
#define UNKFormat 	 0
#define SCFFormat 	 1
#define ABIFormat 	 2
#define MD1Format 	 3
#define MD2Format 	 4

#ifndef FLOAT
#define FLOAT	double
#endif

#ifndef NO_COMPRESS
#define NO_COMPRESS	0
#endif

#ifndef GZ_COMPRESS
#define GZ_COMPRESS	1
#endif

#ifndef Z_COMPRESS
#define Z_COMPRESS	2
#endif

#ifndef BZ_COMPRESS
#define BZ_COMPRESS	3
#endif

#ifndef PATHSEP
#define PATHSEP '/'
#endif

#ifdef PHRED_PATH_MAX
#undef PHRED_PATH_MAX
#endif
#define PHRED_PATH_MAX	4096

#define PHRED_MAX_STRING_LEN    1024

#ifndef int1
#define int1	char
#endif

#ifndef int2
#define int2	short
#endif

#ifndef int4
#define int4	int
#endif

#define MAX_TRACE_VAL	65280.0

typedef unsigned int1	uint1;
typedef unsigned int2	uint2;
typedef unsigned int4	uint4;

#ifndef _WIN32
#define FILE_STATUS	stat
#define STAT_STRUCT	stat
#else
#define FILE_STATUS	_stat
#define STAT_STRUCT	_stat
#endif

