/*********

File		dbReduce.c
Author		EPW <estebanpw@uma.es>
Description	Cuts out a list of genomes from a database to reduce the comparison

USAGE		<genomes_database>		The .fasta database containing the genomes
			<sequence_hits_histogram>	The binary file containing the sequences histogram of hits produced by dictMgMem

**********/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include "metag_common.h"

//goes [a, b)

#define LINE_BUFFER 5000

int main(int argc, char ** av){
	if(argc != 3) terror("USE: dbReduce <genome_database> <sequence_hits_histogram>");
	FILE * f, * histdb, * dbout;
	

	//Open files
	f = fopen64(av[1], "rt");
	if(f == NULL) terror("Error opening genome DB file");

	histdb = fopen64(av[2], "rb");
	if(histdb == NULL) terror("Could not open histogram of sequences");

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

	/*
	//Only for display of values
	uint64_t i;
	for(i=0;i<seqsRead;i++){
		fprintf(stdout, "%"PRIu64" :->: %"PRIu64"\n", i, seqHist[i]);
	}
	exit(-1);
	*/


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
	return 0;
}


