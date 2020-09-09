/** helpParam.c **/

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
int helpParam( void )
#else
int helpParam()
#endif
{
  printf( "\n" );
  printf( "parameter   argument       default    description\n" );
  printf( "---------   --------       -------    -----------\n" );
  printf( "-if         <filename>     none       read input filenames from file\n" );
  printf( "-id         <dirname>      none       read input files from <dirname>\n" );
  printf( "-zd         <dirname>      path       uncompress program path\n" );
  printf( "-zt         <dirname>      /usr/tmp   uncompress temporary directory\n" );
  printf( "\n" );
  printf( "-st         <type>         fasta      sequence file type (fasta|xbap)\n" );
  printf( "-s          none           nofile     write *.seq sequence file(s)\n" );
  printf( "-s          <filename>     nofile     write <filename> sequence file\n" );
  printf( "-sa         <filename>     none       append sequence files to <filename>\n" );
  printf( "-sd         <dirname>      nofile     write *.seq file(s) to <dirname>\n" );
  printf( "-qt         <type>         fasta      quality file type (fasta|xbap|mix)\n" );
  printf( "-q          none           nofile     write *.qual quality file(s)\n" );
  printf( "-q          <filename>     nofile     write <filename> quality file\n" );
  printf( "-qa         <filename>     none       append quality files to <filename>\n" );
  printf( "-qd         <dirname>      nofile     write *.qual file(s) to <dirname>\n" );
  printf( "-qr         <filename>     nofile     write quality report to <filename>\n" );
  printf( "-p          none           nofile     write *.phd.1 file(s)\n" );
  printf( "-p          <filename>     nofile     write <filename> phd file\n" );
  printf( "-pd         <dirname>      nofile     write *.phd.1 file(s) to <dirname>\n" );
  printf( "-xd         <dirname>      nofile     write *.xml file(s) to <dirname>\n" );
  printf( "-cv         <version>      2          SCF format version (2 or 3)\n" );
  printf( "-cp         <precision>    maxval     SCF data precision in bytes (1 or 2)\n" );
  printf( "-cs         none           no scale   always scale traces in SCF files\n" );
  printf( "-c          none           nofile     write * phred SCF file(s)\n" );
  printf( "-c          <filename>     nofile     write <filename> phred SCF file\n" );
  printf( "-cd         <dirname>      nofile     write * SCF file(s) to <dirname>\n" );
  printf( "-d          none           nofile     write *.poly poly file(s)\n" );
  printf( "-d          <filename>     nofile     write <filename> poly file\n" );
  printf( "-dd         <dirname>      nofile     write *.poly file(s) to <dirname>\n" );
  printf( "-raw        <seq name>     NULL       seq name written in output files\n" );
  printf( "-log                       nolog      write phred.log file\n" );
  printf( "\n" );
#ifdef HAVE_GRAPHICS
  printf( "-nocall     none           call       disable basecalling; view ABI\n" );
#else
  printf( "-nocall     none           call       disable basecalling\n" );
#endif
  printf( "-trim       <enzyme seq>   notrim     enable auto trim\n" );
  printf( "-trim_alt   <enzyme seq>   notrim     enable alternate auto trim\n" );
  printf( "-trim_cutoff <n>           0.05       trim_alt error probability\n" );
  printf( "-trim_fasta none           none       trim FASTA bases and qual. values\n" );
  printf( "-trim_scf   none           none       trim SCF bases and qual. values\n" );
  printf( "-trim_phd   none           none       trim base call data in phd files\n" );
  printf( "-trim_out   none           none       trim data in most output files\n" );
  printf( "-nonorm     none           normalize  disable trace normalization\n" );
  printf( "-nosplit    none           none       no compressed peak splitting\n" );
  printf( "-nocmpqv    none           none       no compressed peak quality values\n" );
  printf( "-ceilqv     <ceiling qv>   none       quality value ceiling value\n" );
  printf( "-beg_pred   <point>        none       set peak prediction start point\n" );
  printf( "\n" ) ;
#ifdef HAVE_GRAPHICS
  printf( "-view       none           noview     invoke viewer\n" );
  printf( "-mag        <n>                       display magnification\n" );
  printf( "-bottom     none           top        display sequence complement\n" );
  printf( "-loc        <sub seq>      NULL       locate and display subsequence\n" );
  printf( "-bn         <n>            0          base number to center in display\n" );
  printf( "\n" );
#endif
  printf( "-v          <n>            none       verbose operation <n> = 1 to 63\n" );
  printf( "-tags       none           not tags   label common messages with tags\n" );
  printf( "-V          none           none       show version\n" );
  printf( "-help       none           none       help\n" );
  printf( "-h          none           none       help\n" );
  printf( "-doc        none           none       show phred documentation\n" );
  printf( "\n" );
  printf( "For the warning messages `unable to identify chemistry and dye' and\n" );
  printf( "`unknown chemistry (...) in chromat ...' please read the phred\n" );
  printf( "documentation using the command `phred -doc'.\n" );
  printf( "\n" );

  return( 0 );
}
