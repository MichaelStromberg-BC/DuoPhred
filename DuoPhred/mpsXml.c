/** mpsXml.c **/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "phred.h"

int writeXmlData(int numBases, int *quality, int* qualityUncalled, BaseQual* baseQual, PhredData *phredData) {

	int i;
	int numBytes;
	int numWritten;
	
	FILE *fp;
	
	char *fileRootName;
	char filename[PHRED_PATH_MAX];
	char sequenceName[PHRED_PATH_MAX];
    char buffer[32];
	
	char *bases;
	char *qualities;
    char *positions;
	
	Option *option;

	struct tm *localTimeTm;
	time_t localTimeT;
	char date[64];

	// try to create a new filename
	option = getOption();

	// strip path information
	if((fileRootName = strrchr(option->readFileName, PATHSEP)) != NULL) {
		
		fileRootName++;
	
	} else {
 
		printf("ERROR: Cannot create a new filename in writeXmlData.\n");
		exit(1);
	}

	// construct a sequence filename
	sprintf(sequenceName,"%s",fileRootName);
	fileRootName = strrchr(sequenceName,'.');
	*fileRootName = 0x0;
	
	// construct a new filename
	sprintf(filename,"%s%c%s.xml",option->xmlDataDirName,PATHSEP,sequenceName);

	//
	// construct a new sequence
	//
	numBytes = numBases * sizeof(char) + 1;

	bases = (char*)malloc(numBytes);
	if(bases == NULL) {
		fprintf(stderr, "writeXmlData: error: unable to allocate memory for 'bases'\n" );
		return(ERROR);
	}
  
	for(i=0; i<numBases; i++) {
		bases[i] = "ACGTN"[baseQual[i].nuc];
	}

	bases[numBases] = 0x0;

	//
	// construct a new base quality sequence
	//
	numBytes = (numBases * 3 * sizeof(char)) + 1; // figure on base qualities having two characters

	qualities = (char*)malloc(numBytes);
	if(qualities == NULL) {
		fprintf(stderr, "writeXmlData: error: unable to allocate memory for 'qualities'\n" );
		return(ERROR);
	}

	numWritten = 0;
	for(i=0; i<numBases; i++) {
		numWritten += sprintf(qualities+numWritten,"%d ",quality[i]);
	}

	qualities[numWritten-1] = 0x0;

    //
    // construct a new positions sequence
    //

    // perform a first pass to see how much we need to allocate
    numWritten = 0;
	for(i=0; i<numBases; i++) numWritten += sprintf(buffer,"%d ", baseQual[i].loc);
	
    // allocate the required amount of memory
    numBytes = (numWritten * sizeof(char)) + 1; // figure on base qualities having two characters
	positions = (char*)malloc(numBytes);

    if(positions == NULL) {
		fprintf(stderr, "writeXmlData: error: unable to allocate memory for 'positions'\n" );
		return(ERROR);
	}

    // write our called sequence positions
	numWritten = 0;
	for(i=0; i<numBases; i++) numWritten += sprintf(positions + numWritten, "%d ",baseQual[i].loc);
	positions[numWritten-1] = 0x0;
    
	//
	// write the bases and qualities to the XML file
	//

	// open output file
	if((fp = fopen(filename, "w+")) == NULL) {
		fprintf( stderr, "writeXmlData: error: unable to open file %s\n", filename );
		return(ERROR);
	}

	// construct an XML date
	// TODO: I will need to update the timezone functions
	time(&localTimeT);
	localTimeTm = localtime(&localTimeT);
	sprintf(date,"%4d-%02d-%02dT%02d:%02d:%02d-05:00",localTimeTm->tm_year+1900,localTimeTm->tm_mon+1,localTimeTm->tm_mday,
		localTimeTm->tm_hour,localTimeTm->tm_min,localTimeTm->tm_sec);

	// print the XML header
	fprintf(fp,"<?xml version=\"1.0\" standalone=\"yes\"?>\n");
	fprintf(fp,"<SequenceArchive xmlns=\"http://bayes.bc.edu/SequenceArchive.xsd\">\n");
	
    
    
    fprintf(fp,"  <Sequence Name=\"%s\" Date=\"%s\" ChromatogramFilename=\"%s\" Chemistry=\"",sequenceName, date, getFileName(option->readFileName));

    // write the chemistry type
    if(phredData->chemType == PRIMER_CHEM) fprintf(fp, "primer\" Dye=\""); 
        else if(phredData->chemType == TERMINATOR_CHEM) fprintf(fp, "terminator\" Dye=\""); 
        else fprintf(fp, "unknown\" Dye=\"");

    // write the dye type
    if(phredData->dyeType == RHODAMINE_DYE) fprintf(fp, "rhodamine\" ");
        else if(phredData->dyeType == BIG_DYE_DYE) fprintf(fp, "big\" ");
        else if(phredData->dyeType == ENERGY_TRANSFER_DYE) fprintf(fp, "energy transfer\" ");
        else if(phredData->dyeType == D_RHODAMINE_DYE) fprintf(fp, "d-rhodamine\" ");
        else if(phredData->dyeType == BODIPY_DYE) fprintf(fp, "bodipy\" ");
        else fprintf(fp, "unknown\" ");

    // write the basecaller type
    fprintf(fp,"Basecaller=\"phred %s\">\n", getVersion());

	// write the called bases
	fprintf(fp,"    <Bases>%s</Bases>\n", bases);

    // write the called base qualities
	fprintf(fp,"    <Qualities>%s</Qualities>\n", qualities);

    // write the called base positions
    fprintf(fp,"    <Positions>%s</Positions>\n", positions);
	
	// write the diploid base calls
	fprintf(fp,"    <DiploidBases>\n");
	
	for(i=0; i<numBases; i++) {

		// skip this base if we don't have any uncalled base information
		if(baseQual[i].nucUncalled < 0) continue;
		if(baseQual[i].locUncalled < 0) continue;

		if(qualityUncalled[i] > quality[i]) {
			fprintf(fp,"      <Base index=\"%d\" quality=\"%d\" calledBase=\"%c\" uncalledBase=\"%c\" uncalledPosition=\"%d\" />\n",i,qualityUncalled[i],"ACGTN"[baseQual[i].nuc],"ACGTN"[baseQual[i].nucUncalled],baseQual[i].locUncalled);
		}
	}

	fprintf(fp,"    </DiploidBases>\n");
	fprintf(fp,"  </Sequence>\n");
	fprintf(fp,"</SequenceArchive>\n");

	// close the file
	fclose(fp);

	// delete qualities
	if(qualities) free(qualities);

	// delete bases
	if(bases) free(bases);

	return(OK);
}

