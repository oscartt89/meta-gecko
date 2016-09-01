/*********

File		dbReduce.c
Author		EPW <estebanpw@uma.es>
Description	Cuts out a list of genomes from a database to reduce the comparison

USAGE		<genomes_database>		The .fasta database containing the genomes
			<sequence_hits_histogram>	The binary file containing the sequences histogram of hits produced by dictMgMem
			<working_mode> <n>	Can be [1,2,3...] See below.

WORKING MODES
1	If it has at least n hits [DEFAULT]
2	Select only those genomes with at least n% of hits from the total hits achieved by the most-hitted genome


**********/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include "metag_common.h"
#include "common.h"


//goes [a, b)

#define LINE_BUFFER 5000

int main(int argc, char ** av){
	if(argc != 5) terror("USE: dbReduce <genome_database> <sequence_hits_histogram> <working_mode> <n>");
	FILE * f, * histdb, * dbout;
	

	//Open files
	f = fopen64(av[1], "rt");
	if(f == NULL) terror("Error opening genome DB file");

	histdb = fopen64(av[2], "rb");
	if(histdb == NULL) terror("Could not open histogram of sequences");

	// Working mode selection
	int workingMode = atoi(av[3]);
	uint64_t workingMode_value = asciiToUint64(av[4]);

	//Reading and writing variables
	char buffer[LINE_BUFFER];
	char name_dbout[LINE_BUFFER];
	name_dbout[0] = '\0';



	strcat(name_dbout, av[1]);
	strcat(name_dbout, ".presec");
	fprintf(stdout, "[INFO] Opening output file :%s\n", name_dbout);
	dbout = fopen64(name_dbout, "wt");
	
	//Get number of sequences
    fseeko64(histdb, 0L, SEEK_END);
    uint64_t totalSize = ftello64(histdb);
    fseeko64(histdb, 0L, SEEK_SET);


    //Allocate sequence histogram
    uint64_t nSeqs = totalSize/sizeof(uint64_t);
    fprintf(stdout, "[INFO] Loading %"PRIu64" sequences\n", nSeqs);
	uint64_t * seqHist = (uint64_t *) malloc(nSeqs * sizeof(uint64_t));
	if(seqHist == NULL) terror("Could not allocate memory for sequences histogram");

	//Load sequence histogram into memory
	uint64_t seqsRead = fread(seqHist, sizeof(uint64_t), nSeqs, histdb);
	if(seqsRead != nSeqs) terror("Different amount of read sequences");

	//Create mask to tell which genomes are to be used
	unsigned char * mask = (unsigned char *) malloc(seqsRead*sizeof(unsigned char));
	if(mask == NULL) terror("Could not allocate mask");

	//Fill mask depending on working mode
	uint64_t i;
	switch(workingMode){
		case 2: {
			uint64_t maxFound = 0;
			for(i=0;i<seqsRead;i++){
				if(seqHist[i] > maxFound) maxFound = seqHist[i];
			}
			for(i=0;i<seqsRead;i++){
				mask[i] = (((long double)seqHist[i]/maxFound > (long double) workingMode_value)) ? (1) : (0);
			}
		}
		break;
		//Default is also case 1
		default: {
			for(i=0;i<seqsRead;i++){
				mask[i] = ((seqHist[i] >= workingMode_value)) ? (1) : (0);
			}
		}
	}


	
	//Only for display of values
	for(i=0;i<seqsRead;i++){
		fprintf(stdout, "%"PRIu64" :->: %"PRIu64" [MASK:%d]\n", i, seqHist[i], (int) mask[i]);
	}
	
	


	//Read database
	buffer[0] = 'N'; //Not a '>'
	int64_t genomeCounter=-1;
	while(!feof(f)){
		
		if(buffer[0]=='>'){ //If there is a sequence
			genomeCounter++;

			//If the genome satisfies de condition
			if(seqHist[genomeCounter] > 0){ // TODO put condition here

				fprintf(dbout, "%s", buffer); //Skip first '>'
				fgets(buffer, LINE_BUFFER, f);
				while(buffer[0] != '>' && !feof(f)){
					fprintf(dbout, "%s", buffer);
					fgets(buffer, LINE_BUFFER, f);
				}
			}else{
				fgets(buffer, LINE_BUFFER, f); //Skip the '>' if we dont want it				
			}
		}

		if(buffer[0] != '>') fgets(buffer, LINE_BUFFER, f);
	}
	


	fclose(f);
	fclose(histdb);
	fclose(dbout);

	free(seqHist);
	free(mask);
	return 0;
}


