/** phred.h **/


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


#include <sys/types.h>
#include "typeDef.h"
#include "chromatData.h"
#include "phredData.h"
#include "peak_defs.h"
#include "pstring.h"

#define OK	0
#define ERROR   -1


#define UNKNOWN_CHEM		0
#define PRIMER_CHEM		1
#define TERMINATOR_CHEM		2

#define UNKNOWN_DYE		0
#define RHODAMINE_DYE		1
#define BIG_DYE_DYE		2
#define D_RHODAMINE_DYE		3
#define ENERGY_TRANSFER_DYE	4
#define BODIPY_DYE		5

/*
** Default probability cutoff and minimum sequence
** lenght values used in trimPhred.c and trimAlt.c
*/
#define TRIM_MIN_PRB_VAL     0.05
#define TRIM_MIN_SEQ_LEN     20

/*
** When adding machines, be sure to add the
** machine to the `source' set in chromat2phred.h.
*/
#define UNKNOWN_MACHINE		0
#define ABI_373_377		1
#define MOLDYN_MEGABACE		2
#define	ABI_3700		3
#define LI_COR_4000		4

typedef struct
{
  char *name;
  int   chemType;
  int   dyeType;
  int   machineType;
} ChemistryList;

typedef struct
{
  int numChem;
  int successfulRead;
  ChemistryList *chemList;
} ParFileData;

typedef struct
{
  /*
  ** edit: noedit=0 edit=1
  ** edNum: edit version number
  */
  int edit;
  int editNumber;

  /*
  ** writeSeq: nosequencefile=0  sequencefile=1
  ** sOption: not -s option=0  -s option=1
  ** saOption: not -sa option=0  -sa option=1
  ** sdOption: not -sd option=0  -sd option=1
  ** writeSeqName: *.seq=NULL  <filename>=<string>
  ** seqDirName: seqDirName=<directorname>
  */
  int writeSeq;
  int sOption;
  int saOption;
  int sdOption;
  char *writeSeqName;
  FILE *seqFP;
  char *seqDirName;

  
  /*
  ** writeQual: noqualityfile=0  qualityfile=1
  ** qOption: not -q option=0  -q option=1
  ** qaOption: not -qa option=0  -qa_option=1
  ** qdOption: not -qd option=0  -qd_option=1
  ** writeQualName: *.qual=NULL  <filename>=<string>
  ** qualDirName: qualDirName=<directoryname>
  */
  int writeQual;
  int qOption;
  int qaOption;
  int qdOption;
  char *writeQualName;
  FILE *qualFP;
  char *qualDirName;


  /*
  ** writePhd: nophdfile=0  phd file=1
  ** pOption: not -p option=0  -p option=1
  ** pdOption: not -pd option=0  -pd option=1
  ** writePhdName: *.phd=NULL  <filename>=<string>
  ** phdDirName: phdDirName=<directoryname>
  */
  int writePhd;
  int pOption;
  int pdOption;
  char *writePhdName;
  char *phdDirName;

  /*
  ** writeScf: nofile=0  file=1
  ** writeScfName: *.scf=NULL  <filename>=<string>
  */
  int writeScf;
  char *writeScfName;

  /*
  ** writePolyData: nofile=0  file=1
  ** dOption: not -d option=0  -d option=1
  ** ddOption: not -dd option=0   -dd option=1
  ** writePolyDataName: *.bd=NULL  <filename>=<string>
  ** polyDataDirName: polyDataDirName=<directoryname>
  */
  int writePolyData;
  int dOption;
  int ddOption;
  char *writePolyDataName;
  char *polyDataDirName;

  /*
  ** writeXmlData: nofile=0  file=1
  ** dOption: not -x option=0  -d option=1
  ** xdOption: not -xd option=0   -xd option=1
  ** writeXmlDataName: *.xml=NULL  <filename>=<string>
  ** xmlDataDirName: xmlDataDirName=<directoryname>
  */
  int writeXmlData;
  int xOption;
  int xdOption;
  char *writeXmlDataName;
  char *xmlDataDirName;
    
  /*
  ** bottom
  */
  int bottom;

  /*
  ** trim: 0=notrim 1=trim 2=trim_alt
  */
  int trim;
  char *enzName;

  /*
  ** search and display sequence initially
  */
  int subSeq;
  char *subSequence;

  /*
  ** display magnification (0=default)
  */
  int magnification;

  /*
  ** base number to display initially
  */
  int baseNumber;

  /*
  ** sequence name written in header
  */
  int seq;
  char *seqName;

  /*
  ** call: 0=do not call bases 1=call bases
  */
  int call;

  /*
  ** normalize: 1=normalize traces before calling
  **              bases
  **            0=do not normalize traces
  */
  int normalize;

  /*
  ** seqType: 0=fasta 1=xbap
  */
  int seqType;

  /*
  ** qualType: 0=fasta 1=xbap
  */
  int qualType;

  /*
  ** inType: 0=auto 1=scf 2=abi 3=alf 4=plain
  */
  int inType;

  /*
  ** PHRED SCF output directory
  ** (0=no dir prefix  1=prefix SCF file name with scfDirName)
  */
  int scfDir;
  char *scfDirName;

  /*
  ** input files
  */
  int numInFile;
  int curInFile;
  int ifOption;
  int idOption;
  char *readFileName;
  char **inFileName;
  ino_t inInode;

  /*
  ** enable writing log file
  */
  int log;
  FILE *logfp;

  /*
  ** new/old predicted peak method
  */
  int newpred;

  /*
  ** error flag
  */
  int errorFlag;

  /*
  ** enable diagnostic writing diagnostic files
  */
  int diag;

  /*
  ** Quality report.
  */
  int qualReport;
  char *qualReportFileName;

  /*
  ** SCF data precision.
  */
  int scfPrecOption;
  int scfPrecision;

  /*
  ** SCF output file format.
  */
  int scfVersion;

  /*
  ** Force SCF trace scaling (output with -c options)
  */
  int scaleSCFTraceOption;

  /*
  ** Cross_match discrepancy file read.
  */
  int xmOption;
  char *xmFileName;

  /*
  ** File compression executable directory.
  */
  int filCompressExeDirOption;
  char *filCompressExeDir;

  /*
  ** File compression temporary directory.
  */
  int filCompressTmpDirOption;
  char *filCompressTmpDir; 

  /*
  ** Phred parameter file.
  */
  char parFileName[PHRED_PATH_MAX];
  ParFileData parFileData;

  /*
  ** Compressed peak splitting.
  */
  int compressSplitFlag;

  /*
  ** Primer quality value flag (5 parameter).
  */
  int primerQVFlag;

  /*
  ** Quality value ceiling option.
  */
  int qualityValueCeilingOption;
  int qualityValueCeiling;

  /*
  ** Begin peak prediction point option.
  */
  int beginPeakPredictionOption;
  int beginPeakPredictionPoint;

  /*
  ** Tag option.
  */
  int tagOption;

  /*
  ** Verbose option.
  */
  int verboseOption;
  int verboseLevel;

  /*
  ** Trim point error probability for trim_alt.
  */
  int trimSetOption;
  float trimSetValue;

  /*
  ** Trim FASTA sequence and quality values.
  */
  int trimFastaData;

  /*
  ** Trim SCF sequence and quality values.
  */
  int trimSCFData;

  /*
  ** Trim PHD sequence values.
  */
  int trimPHDData;

} Option;


#ifdef ANSI_C
ChromatData *allocChromatData( int numPoint, int numBase );
int freeChromatData( ChromatData *chromatData );
ChromatData *readData( char *filename, int *status );
PhredData *chromat2phred( ChromatData *chromatData, Option *option );
int autoPhred( Option *option );
int readParam( int argc, char **argv, Option *option );
int traceType( char *filename );
int callSeq( char *inFileName, PhredData *phredData, Option *option, int *status );
int helpParam( void );
int trimSeq( PhredData *phredData, Option *option );
int trimSet( PhredData *phredData );
int writeData( PhredData *phredData, char *filename, Option *option );
int fit_main( int tr_length, FLOAT **tr_vals,
              PhredData *phredData, int *status );
int readRC( Option *option );
int checkParam( Option *option );
int setOption( Option *option );
int initParam( Option *option, int argc );
int openLog( void );
int logParam( int argc, char **argv );
int writeLog( char *string );
int closeLog( void );
int viewPhred( int argc, char **argv, Option *option );
int freeParam( Option *option );
int freeInName( void );
int editOpen( int argc, char **argv, Option *option );
int writeSeq( char *filename, char *seqName, Option *option, PhredData *phredData );
int writeQual( char *filename, char *seqName, Option *option, PhredData *phredData );
int writeSCF( char *filename, int basid, int prec, PhredData *phredData );
int writeHeader( PhredData *phredData, char *seqName, Option *option, FILE *fp, int type );
int setRange( PhredData *phredData );
int readType( char *fn );
float int2float( int in );
int stringMatch( char *seq1, int n1, char *seq2, int n2, int nmiss,
                 int *indices );
char *ourMalloc( int nbytes );
int ourFree( char *ptr );
int editXOpen( int argc, char **argv, Option *option );
Option *getOption( void );
int writePhd( char *filename, char *seqName, Option *option, PhredData *phredData );
int logVersion( void );
char *getVersion( void );
int getMaxVal( FLOAT *maxTrace );
int setMaxVal( int tr_length, FLOAT **tr_vals );
char *getTime( void );
char *getDirName( char *fullPathName );
char *getFileName( char *fullPathName );
int makeTotVals( FLOAT **tr_vals, FLOAT *tot_vals, int tr_length );
int initReport( void );
int makeReport( PhredData *phredData );
int writeReport( Option *option );
BaseQual *initBaseQual( Peak *first_peak, int *numBase );
int maxDownBaseQual( int tr_length, FLOAT **tr_vals, FLOAT *tot_vals, int numBase, BaseQual *baseQual );
LocPeak *evaluateTrace( int tr_length, FLOAT **tr_vals, FLOAT *tot_vals, int block_size,
                       int *numLocPeak, int *begGoodTrace, int *endGoodTrace );
int max_from_fft( FLOAT *vec, int length, FLOAT *period, FLOAT *value, FLOAT start_check, FLOAT end_check );
int findStartPred( int tr_length, int npk, LocPeak *locPeak, int block_size, int *tr_qual );
int makePeakSignal( int tr_length, int npk, LocPeak *locPeak, FLOAT *waveForm );
int init_peak( Peak *peak, FLOAT pred_location, FLOAT pred_period,
               FLOAT total_signal, FLOAT proportion_fitted );
FLOAT nearest_peak( FLOAT point, FLOAT offset, FLOAT period );
int find_best( FLOAT *tot_vals, int length, FLOAT start_period,
               FLOAT min_period, FLOAT max_period, FLOAT period_incr,
               FLOAT *add_max_value, FLOAT *add_best_period, FLOAT *add_best_offset );
FLOAT weight_vec( FLOAT *temp_vec, FLOAT *vec, int length, int site );
int free_cos_sin( void );
FLOAT alignment_score( Observed_peak *obs_peak, int m,
                       int k, Peak *peak, FLOAT shift );
FLOAT measure_shift( FLOAT location, Peak *peak );
int normalize_traces( FILE *fp, FLOAT **tr_vals, FLOAT *tot_vals,
                      int tr_length, FLOAT *scale_traces );
Peak *make_predicted_peaks( FILE *fp, FLOAT **tr_vals, int tr_length,
                            FLOAT *tot_vals, int numLocPeak, LocPeak *locPeak,
                            int beginPeakPredictionOption,
                            int beginPeakPredictionPoint );
Peak *alloc_predicted_peak( void );
Observed_peak *find_observed_peaks( FLOAT **tr_vals, int tr_length );
int fit_peaks( int tr_length, FLOAT **tr_vals, Peak *first_peak,
               Observed_peak *first_obs_peak );
int capture_missed_peaks1( int tr_length, Peak *first_peak,
                          Observed_peak *first_obs_peak, FLOAT **sig_func,
                          int begGoodTrace, int endGoodTrace );
int capture_missed_peaks2( int tr_length, FLOAT **tr_vals, FLOAT *tot_vals,
                          Peak *first_peak, Observed_peak *first_obs_peak );
int capture_missed_peaks3( int tr_length, FLOAT **tr_vals, FLOAT *tot_vals,
                           Peak *first_peak, Observed_peak *first_obs_peak,
                           int npk, LocPeak *locPeak );
Peak *fit_sine( FLOAT **tr_vals, FLOAT *tot_vals, int tr_length,
                FLOAT *scale_traces, Observed_peak **pfirst_obs_peak,
                int compressSplitFlag, int chemType, int *status );
int setQual( int tr_length, FLOAT **tr_vals, FLOAT *tot_vals,
             int numBase, Peak *first_peak, int *quality, int *qualityUncalled,
             PhredData *phredData, int qvCeilingOption, int qvCeiling );
int fit_compress( Peak *first_peak, int tr_length, FLOAT **tr_vals, FLOAT *tot_vals );
int loadPolyData( Peak *first_peak, FLOAT *scale_traces, FLOAT **tr_vals,
                  PolyData *polyData );
int freePolyData( PolyData *polyData );
int free_obs_peak( Observed_peak *first_obs_peak );
int free_predicted_peak( void );
int freePhredData( PhredData *phredData );
Observed_peak *make_observed_peak( int nuc );
int freeFOF( void );
int freeDIR( void );
int initFOF( void );
int initDIR( void );
Cos_sin_array *make_cos_sin( double period, int length );
double inner( FLOAT *vec, int length, double period, double *sum_sin );
int readFOF( int argc, char *filename, Option *option );
int readDIR( int argc, char *dirname, Option *option );
int spcRatioBaseQual( int numBase, BaseQual *baseQual );
int maxRatioBaseQual( int tr_length, FLOAT **tr_vals, FLOAT **otr_vals,
                      int un_norm, int numBase, BaseQual *baseQual );
FLOAT get_area( FLOAT *trace, int startp, int endp );
FLOAT max_area( FLOAT *tx, FLOAT *ty, FLOAT *tz, int stp, int endp );
int writePolyData( char *filename, char *seqname, int numPoint, PolyData *polyData );
int viewOpen( int argc, char **argv, Option *option );
int writeSCF2( char *fn, int prec, PhredData *phredData );
int writeSCF3( char *fn, int prec, PhredData *phredData );
int freeChemList( int numChem, ChemistryList *chemList );

int trimAlt( PhredData *phredData, char *enzName );
int trimPhred( PhredData *phredData );
int compressMotif( int numBase, char *seq, int *compress );
int resPar7( int tr_length, FLOAT **tr_vals, FLOAT *tot_vals,
             int numBase, BaseQual *baseQual );
int readEnvVar( Option *option );
int readParamFile( Option *option );
int delFile( char *filename );
int compressParameter( int numBase, BaseQual *baseQual );
int testCompress( char *ifnm, char *ofnm );
int uncompressFile( char *ifnm, char *ofnm, Option *option, int type );
int reportFileStatus( Option *option, char *inFileName, int statusCode );
int isESD( char *filename );
int compareFiles( char *filename, Option *option );
int showDoc( void );
int findTraceExtrema( ChromatData *chromatData );
FLOAT getMaxRatioN(int loc, int windowLen, FLOAT **tr_vals, BaseQual *baseQual);
FLOAT getMaxRatioNU(int loc, int windowLen, FLOAT **tr_vals, BaseQual *baseQual);
int compare(void *context, const void *v1, const void *v2);
int writeXmlData(int numBases, int *quality, int* qualityUncalled, BaseQual* baseQual, PhredData *phredData);

#else
ChromatData *readABI();
ChromatData *readSCF();
ChromatData *allocChromatData();
int freeChromatData();
ChromatData *readData();
PhredData *chromat2phred();
int autoPhred();
int readParam();
int traceType();
int callSeq();
int helpParam();
int trimSeq();
int trimSet();
int writeData();
int fit_main();
int readRC();
int checkParam();
int setOption();
int initParam();
int openLog();
int logParam();
int writeLog();
int closeLog();
int viewPhred();
int freeParam();
int freeInName();
int editOpen();
int writeSeq();
int writeQual();
int writeSCF();
int writeHeader();
int setRange();
int readType();
float int2float();
int stringMatch();
char *ourMalloc();
int ourFree();
int editXOpen();
Option *getOption();
int writePhd();
int logVersion();
char *getVersion();
int getMaxVal();
int setMaxVal();
char *getTime();
char *getDirName();
char *getFileName();
int makeTotVals();
int initReport();
int makeReport();
int writeReport();
BaseQual *initBaseQual();
int maxDownBaseQual();
LocPeak *evaluateTrace();
int max_from_fft();
int findStartPred();
int makePeakSignal();
int init_peak();
FLOAT nearest_peak();
int find_best();
FLOAT weight_vec();
int free_cos_sin();
FLOAT alignment_score();
FLOAT measure_shift();
Peak *make_predicted_peaks();
Peak *alloc_predicted_peak();
Observed_peak *find_observed_peaks();
int fit_peaks();
int capture_missed_peaks1();
int capture_missed_peaks2();
int capture_missed_peaks3();
Peak *fit_sine();
int setQual();
int fit_compress();
int loadPolyData();
int freePolyData();
int free_obs_peak();
int free_predicted_peak();
int freePhredData();
Observed_peak *make_observed_peak();
int freeFOF();
int freeDIR();
int initFOF();
int initDIR();
Cos_sin_array *make_cos_sin();
double inner();
int readFOF();
int readDIR();
int spcRatioBaseQual();
int maxRatioBaseQual();
FLOAT get_area();
FLOAT max_area();
int writePolyData();
int viewOpen();
int writeSCF2();
int writeSCF3();
int freeChemList();

int trimAlt();
int trimPhred();
int compressMotif();
int resPar7();
int readEnvVar();
int readParamFile();
int delFile();
int compressParameter();
int testCompress();
int uncompressFile();
int reportFileStatus();

int isESD();
int compareFiles();
int showDoc();
int findTraceExtrema();

#endif

